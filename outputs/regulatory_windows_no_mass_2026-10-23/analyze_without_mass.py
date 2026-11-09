"""Refit the six-window regulatory peak-count analysis without body mass.

Reuse validated selected-locus counts without modifying the measured inputs.
Run with /usr/bin/python3 (NumPy and SciPy are required).
"""
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, linregress, pearsonr, t

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
SOURCE = ROOT / "outputs/regulatory_windows_2026-10-22"
WINDOWS = [20, 50, 100, 200, 500, 1000]
FIELDS = ["window_kb", "gene", "n_species", "r", "p_value", "BH_FDR_window",
          "BH_FDR_all_windows", "MLS_coefficient", "intercept", "df", "species",
          "direction", "mass_adjusted_r", "mass_adjusted_FDR"]


def read(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path, rows, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    source_manifest = json.loads((SOURCE / "manifest.json").read_text())
    for source in source_manifest["source_files"]:
        path = Path(source["path"])
        stat = path.stat() if path.exists() else None
        actual = (stat.st_size, stat.st_mtime_ns) if stat else (None, None)
        assert actual == (source["size"], source["mtime_ns"]), f"Changed input: {path}"
    database = Path(source_manifest["database"])
    assert database.stat().st_mtime_ns == source_manifest["database_mtime_ns"]
    assert hashlib.sha256((database.parent / "species_mls.json").read_bytes()).hexdigest() == source_manifest["traits_sha256"]
    species = read(SOURCE / "species.tsv")
    # There are no extra species eligible only after removing the mass requirement.
    assert not any(r["analysis_status"] == "Body mass unavailable" for r in species)
    traits = {r["species"]: float(r["MLS_years"]) for r in species if r["analysis_status"] == "Eligible"}
    prior_summary = {int(r["window_kb"]): r for r in read(SOURCE / "window_summary.tsv")}
    results, summaries, validation = {}, [], []
    input_hashes = {}

    for kb in WINDOWS:
        measurements = defaultdict(dict)
        for sid in sorted(traits):
            path = SOURCE / "selected_loci" / f"{sid}_{kb}kb.tsv.gz"
            input_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
            records = read(path)
            assert len({r["gene"] for r in records}) == len(records)
            for record in records:
                assert record["species"] == sid
                eligible = float(record["available_fraction"]) >= 0.5
                assert eligible == (record["region_eligible"] == "True")
                count = float(record["median"])
                if eligible and math.isfinite(count) and count > 0:
                    measurements[record["gene"]][sid] = count

        old = {r["gene"]: r for r in read(SOURCE / f"tss_{kb}kb/all_tested.tsv.gz")}
        rows, vectors = [], {}
        for gene, values in sorted(measurements.items()):
            ids = sorted(values)
            n = len(ids)
            if n < 5:
                continue
            x = np.log10([traits[s] for s in ids])
            y = np.log10([values[s] for s in ids])
            if np.ptp(x) == 0 or np.ptp(y) == 0:
                continue
            xc, yc = x-x.mean(), y-y.mean()
            xx, yy, xy = float(xc@xc), float(yc@yc), float(xc@yc)
            r = float(np.clip(xy/math.sqrt(xx*yy), -1, 1))
            df = n-2
            p = float(2*t.sf(abs(r)*math.sqrt(df/max(1-r*r, 1e-300)), df))
            prior = old.get(gene)
            if prior:
                assert prior["species"] == ",".join(ids)
                assert int(prior["n_species"]) == n
            row = dict(window_kb=kb, gene=gene, n_species=n, r=r, p_value=p,
                       MLS_coefficient=xy/xx, intercept=float(y.mean()-(xy/xx)*x.mean()),
                       df=df, species=",".join(ids), direction="Positive" if r > 0 else "Negative" if r < 0 else "Zero",
                       mass_adjusted_r=float(prior["partial_r"]) if prior else None,
                       mass_adjusted_FDR=float(prior["BH_FDR_window"]) if prior else None)
            rows.append(row)
            vectors[gene] = (x, y)
        assert set(old) <= {r["gene"] for r in rows}, "Previously testable genes were lost"

        q = false_discovery_control([r["p_value"] for r in rows], method="bh")
        # Independent BH implementation validates the correction family and values.
        order = np.argsort([r["p_value"] for r in rows])
        sorted_p = np.array([rows[i]["p_value"] for i in order])
        q_manual = np.minimum.accumulate((sorted_p*len(rows)/np.arange(1, len(rows)+1))[::-1])[::-1]
        assert np.allclose(q[order], np.minimum(q_manual, 1), rtol=1e-12, atol=1e-15)
        for row, value in zip(rows, q):
            row["BH_FDR_window"] = float(value)
        rows.sort(key=lambda r: (r["p_value"], r["gene"]))

        # Independently check all discoveries plus a deterministic nonsignificant sample.
        chosen = {r["gene"] for r in rows if r["BH_FDR_window"] < 0.05}
        chosen.update(rows[i]["gene"] for i in np.linspace(0, len(rows)-1, min(100, len(rows)), dtype=int))
        chosen.update(g for g in ["SIRT6", "FBXO22", "FBXO25"] if g in vectors)
        lookup = {r["gene"]: r for r in rows}
        max_r_diff, max_p_diff = 0., 0.
        for gene in sorted(chosen):
            x, y = vectors[gene]
            reference_r, reference_p = pearsonr(x, y)
            fit = linregress(x, y)
            row = lookup[gene]
            max_r_diff = max(max_r_diff, abs(row["r"]-reference_r))
            max_p_diff = max(max_p_diff, abs(row["p_value"]-reference_p), abs(row["p_value"]-fit.pvalue))
            assert math.isclose(row["r"], reference_r, abs_tol=1e-12)
            assert math.isclose(row["p_value"], reference_p, rel_tol=1e-9, abs_tol=1e-12)
            assert math.isclose(row["p_value"], fit.pvalue, rel_tol=1e-9, abs_tol=1e-12)
            assert math.isclose(row["MLS_coefficient"], fit.slope, rel_tol=1e-10, abs_tol=1e-12)
        validation.append(dict(window_kb=kb, independent_fits=len(chosen), maximum_r_difference=float(max_r_diff),
                               maximum_p_difference=float(max_p_diff), same_tested_genes=set(old)==set(lookup),
                               added_testable_genes=sorted(set(lookup)-set(old))))
        results[kb] = rows
        sig = [r for r in rows if r["BH_FDR_window"] < 0.05]
        prior_sig = {g for g, r in old.items() if float(r["BH_FDR_window"]) < .05}
        current_sig = {r["gene"] for r in sig}
        summary = dict(window_kb=kb, full_width_kb=2*kb, tested=len(rows),
                       FDR_positive=sum(r["r"] > 0 for r in sig), FDR_negative=sum(r["r"] < 0 for r in sig), FDR_total=len(sig),
                       nominal_positive=sum(r["r"] > 0 and r["p_value"] < .05 for r in rows),
                       nominal_negative=sum(r["r"] < 0 and r["p_value"] < .05 for r in rows),
                       mass_adjusted_FDR_positive=int(prior_summary[kb]["FDR_positive"]),
                       mass_adjusted_FDR_negative=int(prior_summary[kb]["FDR_negative"]),
                       shared_FDR_genes=len(current_sig & prior_sig), new_FDR_genes=len(current_sig-prior_sig),
                       lost_FDR_genes=len(prior_sig-current_sig), minimum_FDR=min(r["BH_FDR_window"] for r in rows))
        summaries.append(summary)
        print(json.dumps(summary), flush=True)

    all_rows = [r for rows in results.values() for r in rows]
    for row, value in zip(all_rows, false_discovery_control([r["p_value"] for r in all_rows], method="bh")):
        row["BH_FDR_all_windows"] = float(value)
    all_significant = []
    for summary in summaries:
        kb = summary["window_kb"]
        rows = results[kb]
        folder = OUT/f"tss_{kb}kb"
        folder.mkdir(exist_ok=True)
        write(folder/"all_tested.tsv.gz", rows, FIELDS)
        for direction, sign in [("positive", 1), ("negative", -1)]:
            selected = [r for r in rows if r["r"]*sign > 0 and r["BH_FDR_window"] < .05]
            write(folder/f"FDR05_{direction}.tsv", selected, FIELDS)
            (folder/f"FDR05_{direction}_genes.txt").write_text("".join(r["gene"]+"\n" for r in selected))
        summary["joint_FDR_positive"] = sum(r["r"] > 0 and r["BH_FDR_all_windows"] < .05 for r in rows)
        summary["joint_FDR_negative"] = sum(r["r"] < 0 and r["BH_FDR_all_windows"] < .05 for r in rows)
        all_significant.extend(r for r in rows if r["BH_FDR_window"] < .05)
    write(OUT/"window_summary.tsv", summaries, list(summaries[0]))
    write(OUT/"significant_genes_all_windows.tsv", all_significant, FIELDS)
    write(OUT/"species.tsv", species, list(species[0]))
    union = sorted({r["gene"] for r in all_significant})
    (OUT/"significant_genes_union.txt").write_text("".join(g+"\n" for g in union))
    manifest = dict(source_directory=str(SOURCE), model="log10(species median ATAC peak count) ~ log10(MLS years)",
                    covariates=[], measurement="ATAC peak count, not integrated CPM within peaks",
                    windows_flank_kb=WINDOWS, significance="BH FDR < 0.05 separately within each window",
                    additional_significance="BH across all gene-window tests also supplied",
                    species_eligible=len(traits), excluded_species=source_manifest["excluded_species"],
                    lower_bound_MLS_excluded=True, minimum_species=5, minimum_available_fraction=.5,
                    locus_selection=source_manifest["locus_selection"], availability=source_manifest["availability"],
                    source_files_verified=len(source_manifest["source_files"]), annotation_and_traits_verified=True,
                    selected_locus_sha256=input_hashes, validation=validation, total_tests=len(all_rows),
                    significant_gene_window_pairs=len(all_significant), unique_significant_genes=len(union), summary=summaries)
    (OUT/"manifest.json").write_text(json.dumps(manifest, indent=2))
    (OUT/"README.txt").write_text(
        "Regulatory activity versus MLS without body-mass correction — 23 October 2026\n\n"
        "Model: log10(species-median ATAC peak count) ~ log10(MLS in years). No covariates.\n"
        "The 20, 50, 100, 200, 500 and 1,000 kb windows are symmetric flanks around the TSS; their full widths are twice these values.\n"
        "Primary significance: BH FDR < 0.05 separately across all tested genes within each window. Joint BH across six windows is provided separately.\n"
        "Uses the same peak counts, highest-median locus selection, traits and available species as the 22 October mass-adjusted analysis. All 150 peak/bigWig source-file metadata records, annotation database mtime and traits hash were verified unchanged.\n"
        "Human, macaque, bowhead whale and missing/lower-bound MLS are excluded. At least five species per gene. Nonpositive medians are excluded without a pseudocount. Available window span must be >=50%; availability measures contig bounds, not internal sequence gaps or ATAC coverage. Highest-median locus is selected before filtering, with no fallback.\n"
        "The predictor is log10(MLS), not raw MLS. Regulatory activity here means peak count, matching the previous six-window analysis; it does not mean integrated CPM restricted to peaks.\n"
        "No phylogenetic, body-mass, sequencing-depth, batch or global peak-count correction is applied. Species may differ between genes/windows. The original mass-adjusted analysis and browser defaults were not modified.\n\n"
        "window_summary.tsv: counts by window, directions and comparison with mass adjustment.\n"
        "significant_genes_all_windows.tsv: all within-window FDR-significant gene-window associations.\n"
        "tss_*kb/: all tested statistics and separate positive/negative significant gene lists.\n"
        "manifest.json: settings, source hashes and independent SciPy/BH validation.\n")
    print(json.dumps({"unique_significant_genes": len(union), "gene_window_pairs": len(all_significant),
                      "tests": len(all_rows), "validation": validation}), flush=True)


if __name__ == "__main__":
    main()
