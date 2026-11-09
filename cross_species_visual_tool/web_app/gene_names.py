"""Shared, offline gene-symbol/ortholog matching without merging paralogs."""
import gzip
import json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).resolve().parent/'data/gene_names.json.gz'
LABELS = {'exact_symbol':'Exact symbol', 'symbol_alias':'HGNC symbol alias',
          'ortholog':'HCOP ortholog', 'projected_ortholog':'Mouse-derived annotation',
          'ambiguous':'Ambiguous gene mapping', 'not_found':'No supported gene match'}


class GeneNames:
    def __init__(self, data):
        self.data = data
        self.approved = set(data.get('approved', []))
        self.aliases = data.get('aliases', {})
        self.reverse_aliases = defaultdict(set)
        for alias, canonical in self.aliases.items():
            self.reverse_aliases[canonical].add(alias)
        self.orthologs = data.get('orthologs', {})
        self.inputs = defaultdict(set)
        for mapping in self.orthologs.values():
            for canonical, entries in mapping.items():
                for entry in entries:
                    if entry['reciprocal_unique']:
                        self.inputs[entry['symbol'].upper()].add(canonical)
        self.projected = {sid:set(names) for sid,names in data.get('projected_mouse_symbols', {}).items()}

    def canonical(self, query):
        key = str(query).strip().upper()
        if key in self.approved:
            return key
        choices = self.inputs.get(key, set()) | ({self.aliases[key]} if key in self.aliases else set())
        return next(iter(choices)) if len(choices) == 1 else key

    def candidates(self, query):
        canonical = self.canonical(query)
        keys = {canonical} | self.reverse_aliases.get(canonical, set())
        for mapping in self.orthologs.values():
            keys.update(e['symbol'].upper() for e in mapping.get(canonical, []))
        return sorted(keys)

    def resolve(self, query, sid, available, allow_projected=True):
        canonical = self.canonical(query)
        if not isinstance(available, (dict, set)):
            available = set(available)
        native = self.data.get('native_symbols', {}).get(sid, {})

        def result(method, symbol=None, candidates=None, support=None):
            source = ('https://www.genenames.org/tools/hcop/#!/?q='+canonical+'&qtype=symbol&qtax_id=9606&submit=true'
                      if method in ('ortholog','projected_ortholog','ambiguous') else
                      'https://www.genenames.org/tools/multi-symbol-checker/' if method == 'symbol_alias' else None)
            return dict(query=str(query).strip().upper(), canonical=canonical, symbol=symbol,
                        native_symbol=native.get(symbol, symbol), method=method, label=LABELS[method],
                        candidates=sorted(candidates or []), support=support or [], source=source)

        if canonical in available:
            return result('exact_symbol', canonical)
        evidence = []
        for entry in self.orthologs.get(sid, {}).get(canonical, []):
            if entry['symbol'].upper() in available:
                evidence.append((entry, 'ortholog'))
        if allow_projected and sid in self.projected:
            for entry in self.orthologs.get('Mouse', {}).get(canonical, []):
                key = entry['symbol'].upper()
                if key in available and key in self.projected[sid]:
                    evidence.append((entry, 'projected_ortholog'))
        if evidence:
            keys = {e['symbol'].upper() for e,_ in evidence}
            if len(keys) != 1 or any(not e['reciprocal_unique'] for e,_ in evidence):
                return result('ambiguous', candidates=keys)
            entry, method = evidence[0]
            return result(method, entry['symbol'].upper(), support=entry['support'])
        keys = {key for key in self.reverse_aliases.get(canonical, set()) if key in available}
        if len(keys) == 1:
            return result('symbol_alias', next(iter(keys)))
        if keys:
            return result('ambiguous', candidates=keys)
        return result('not_found')


@lru_cache(maxsize=1)
def resolver():
    with gzip.open(DATA, 'rt') as handle:
        return GeneNames(json.load(handle))


def canonical(query):
    return resolver().canonical(query)


def resolve(query, sid, available, allow_projected=True):
    return resolver().resolve(query, sid, available, allow_projected)


def candidates(query):
    return resolver().candidates(query)


def match_records(query, records, species_key='species', symbol_key='symbol_key', allow_projected=True):
    """Select one supported symbol per species. Keep separate loci and samples."""
    grouped = defaultdict(list)
    for row in records:
        grouped[row[species_key]].append(row)
    latin_ids = {latin:sid for sid,latin in resolver().data.get('species_latin', {}).items()}
    selected, matches = [], {}
    for species, rows in grouped.items():
        sid = latin_ids.get(species, species)
        match = resolve(query, sid, {r[symbol_key].upper() for r in rows}, allow_projected)
        matched_rows = [r for r in rows if r[symbol_key].upper() == match['symbol']]
        if matched_rows:
            match['native_symbol'] = matched_rows[0].get('symbol', match['symbol'])
        matches[species] = match
        for row in rows:
            if row[symbol_key].upper() == match['symbol']:
                selected.append(dict(row, gene_match=match))
    return selected, matches
