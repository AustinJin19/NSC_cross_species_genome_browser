#!/usr/bin/env python3
"""Build adult liver SE BEDs from the saved dbSUPER snapshot (Python stdlib only)."""
import bisect
import collections
import csv
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlencode

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "sources"
PROJECT = ROOT.parent.parent
BASE = "https://dbsuper.khanlab.bio/api/"
MALIGNANT = re.compile(
    r"carcinoma|adenoma|cancer|sarcoma|leukemia|lymphoma|melanoma|glioma|"
    r"myeloma|malignan|neoplasm|tumou?r", re.I
)
CONFIG = {
    "human_hg38": {
        "assembly": "GRCh38", "organism": "Homo sapiens",
        "fai": "reference_genomes/selected3/hg38.fa.fai",
        "liver_ids": {"EXH0166", "EXH0282", "EXH0320", "EXH0340", "EXH0368"},
    },
    "mouse_mm10": {
        "assembly": "GRCm38", "organism": "Mus musculus",
        "fai": "reference_genomes/mouse_mm10/mm10.fa.fai",
        "liver_ids": {"EXM0064"},
    },
}


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def catalog(stem):
    pages = [json.loads((SRC / f"{stem}{suffix}.json").read_text())
             for suffix in ("", "_p2", "_p3")]
    rows = [r for page in pages for r in page["results"]]
    assert len(rows) == pages[0]["count"] == len({r["identifier"] for r in rows})
    return {r["identifier"]: r for r in rows}


def chrom_key(chrom):
    suffix = chrom.removeprefix("chr")
    if suffix.isdigit():
        return (0, int(suffix), "")
    return (1, {"X": 0, "Y": 1, "M": 2}.get(suffix, 3), suffix)


def key(row):
    return chrom_key(row["chrom"]), row["start"], row["end"], row["name"]


def read_bed(path):
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track ", "browser ")):
            continue
        c, a, b, name, *rest = line.split("\t")
        assert re.fullmatch(r"EX[HM]\d{4}SE\d+", name), (path, line)
        rows.append(dict(chrom=c, start=int(a), end=int(b), name=name,
                         experiment=name.split("SE")[0],
                         source_score=float(rest[0]) if rest else None))
    assert len(rows) == len({r["name"] for r in rows}), path
    return rows


def write_bed(path, rows):
    # BED4 deliberately omits raw ROSE scores: these exceed BED's 0--1000 range.
    path.write_text("".join(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['name']}\n"
                            for r in sorted(rows, key=key)))


def merge_calls(rows):
    """Union strictly overlapping spans; do not join merely adjacent intervals."""
    result = []
    for row in sorted(rows, key=key):
        if result and row["chrom"] == result[-1]["chrom"] and row["start"] < result[-1]["end"]:
            result[-1]["end"] = max(result[-1]["end"], row["end"])
            result[-1]["calls"].append(row)
        else:
            result.append(dict(chrom=row["chrom"], start=row["start"],
                               end=row["end"], calls=[row]))
    return result


def overlap_index(rows):
    by_chrom = collections.defaultdict(list)
    for r in merge_calls(rows):
        by_chrom[r["chrom"]].append((r["start"], r["end"]))
    return {chrom: ([a for a, b in spans], [b for a, b in spans])
            for chrom, spans in by_chrom.items()}


def overlaps(row, index):
    starts, ends = index.get(row["chrom"], ([], []))
    i = bisect.bisect_right(ends, row["start"])
    return i < len(starts) and starts[i] < row["end"]


def qualifies(exp, sample, cfg):
    return (sample["organism"] == cfg["organism"]
            and any(a["name"] == cfg["assembly"] for a in exp["assemblies"])
            and exp["target"] == "H3K27ac" and not exp["is_perturbed"]
            and sample["sample_class"] == "tissue"
            and sample["life_stage"].lower().startswith("adult")
            and not MALIGNANT.search(" ".join(str(sample.get(k, ""))
                                            for k in ("name", "description", "disease"))))


def is_liver(sample):
    return "liver" in (sample["name"] + " " + sample["tissue"]).lower()


def main():
    samples = catalog("all_samples")
    experiments = catalog("all_experiments")
    summary = {}
    all_datasets = []
    generated = []
    for prefix, cfg in CONFIG.items():
        chosen = {i: e for i, e in experiments.items()
                  if qualifies(e, samples[e["sample"]], cfg)}
        liver_ids = {i for i, e in chosen.items() if is_liver(samples[e["sample"]])}
        assert liver_ids == cfg["liver_ids"], liver_ids
        comparator_ids = set(chosen) - liver_ids
        if prefix == "human_hg38":
            calls = []
            for r in csv.DictReader((SRC / "human_liver_native.csv").open()):
                assert r["Genome"] == cfg["assembly"]
                assert int(r["Size (bp)"]) == int(r["End"]) - int(r["Start"])
                calls.append(dict(chrom=r["Chromosome"], start=int(r["Start"]),
                                  end=int(r["End"]), name=r["Super Enhancer ID"],
                                  experiment=r["Experiment"], source_score=float(r["Score"])))
            native = read_bed(SRC / "human_liver_native.bed")
            coordinates = lambda rr: {(r["name"], r["chrom"], r["start"], r["end"]) for r in rr}
            assert coordinates(native) == coordinates(calls), "Native BED/CSV coordinate mismatch"
            requests = json.loads((SRC / "comparator_download_requests.json").read_text())
            comparator_calls = []
            for req in requests:
                chunk = read_bed(SRC / req["file"])
                reported = json.loads((SRC / req["file"].replace(".bed", ".count.json")).read_text())
                assert len(chunk) == reported["export_rows"], req["file"]
                comparator_calls.extend(chunk)
        else:
            calls = read_bed(SRC / "mouse_liver_native.bed")
            comparator_calls = [row for row in read_bed(SRC / "mouse_tissue_native.bed")
                                if row["experiment"] in comparator_ids]
        assert {r["experiment"] for r in calls} == liver_ids
        assert {r["experiment"] for r in comparator_calls} == comparator_ids
        assert len(comparator_calls) == len({r["name"] for r in comparator_calls})
        counts = collections.Counter(r["experiment"] for r in calls + comparator_calls)

        sizes_path = SRC / f"{prefix}.chrom.sizes"
        if not sizes_path.exists():
            sizes_path.write_text("".join("\t".join(line.split("\t")[:2]) + "\n"
                                         for line in (PROJECT / cfg["fai"]).read_text().splitlines()))
        sizes = {c: int(n) for c, n in (line.split("\t") for line in sizes_path.read_text().splitlines())}
        for r in calls + comparator_calls:
            assert r["chrom"] in sizes and 0 <= r["start"] < r["end"] <= sizes[r["chrom"]], r
        index = overlap_index(comparator_calls)
        merged = merge_calls(calls)
        for n, r in enumerate(merged, 1):
            r["name"] = f"{prefix}_liver_SE_{n:04d}"
            r["liver_datasets"] = sorted({c["experiment"] for c in r["calls"]})
            r["n_liver_datasets"] = len(r["liver_datasets"])
            r["candidate_liver_specific"] = not overlaps(r, index)
        candidates = [r for r in merged if r["candidate_liver_specific"]]
        supported = [r for r in candidates if r["n_liver_datasets"] >= 2]
        outputs = {
            f"{prefix}_liver_super_enhancers.bed": merged,
            f"{prefix}_candidate_liver_specific_super_enhancers.bed": candidates,
        }
        if len(liver_ids) > 1:
            outputs[f"{prefix}_candidate_liver_specific_supported2datasets.bed"] = supported
        for filename, rows in outputs.items():
            write_bed(ROOT / filename, rows)
            generated.append(filename)
        raw_dir = ROOT / "per_dataset"
        raw_dir.mkdir(exist_ok=True)
        for eid in sorted(liver_ids):
            write_bed(raw_dir / f"{prefix}_{eid}.bed", [r for r in calls if r["experiment"] == eid])
        save_json(ROOT / f"{prefix}_region_metadata.json", merged)

        for eid, exp in sorted(chosen.items()):
            sample = samples[exp["sample"]]
            all_datasets.append(dict(
                organism=cfg["organism"], assembly=cfg["assembly"],
                role="liver" if eid in liver_ids else "adult_nonliver_comparator",
                super_enhancer_count=counts[eid], experiment=exp, sample=sample,
                encode_url=f"https://www.encodeproject.org/experiments/{exp['case_acc']}/",
                dbsuper_api_url=BASE + "super-enhancers/?" + urlencode({"experiment": eid, "assembly": cfg["assembly"]}),
            ))
        summary[prefix] = dict(
            assembly=cfg["assembly"], organism=cfg["organism"],
            adult_liver_datasets=len(liver_ids), source_liver_call_count=len(calls),
            merged_liver_regions=len(merged), candidate_liver_specific_regions=len(candidates),
            candidate_supported_by_at_least_two_datasets=len(supported) if len(liver_ids) > 1 else None,
            adult_nonliver_datasets=len(comparator_ids), nonliver_call_count=len(comparator_calls),
            nonliver_tissue_labels=sorted({samples[chosen[e]["sample"]]["tissue"] for e in comparator_ids}),
            liver_call_counts={e: counts[e] for e in sorted(liver_ids)},
            total_merged_liver_bp=sum(r["end"]-r["start"] for r in merged),
            total_candidate_bp=sum(r["end"]-r["start"] for r in candidates),
            chromosome_size_validation="All source liver and selected comparator intervals passed",
        )
    save_json(ROOT / "summary.json", summary)
    save_json(ROOT / "dataset_provenance.json", all_datasets)
    # Source snapshot checksums allow reruns to be distinguished from later database updates.
    source_files = sorted(p for p in SRC.iterdir() if p.is_file() and p.suffix in {".json", ".bed", ".csv", ".sizes"})
    save_json(ROOT / "source_checksums.json", [dict(file=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
              sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in source_files])
    urls = {}
    for stem, endpoint in (("all_samples", "samples"), ("all_experiments", "experiments")):
        for page, suffix in enumerate(("", "_p2", "_p3"), 1):
            urls[f"{stem}{suffix}.json"] = BASE + endpoint + "/?" + urlencode({"format": "json", "page_size": 200, "page": page})
    human_ids = ",".join(sorted(CONFIG["human_hg38"]["liver_ids"]))
    for fmt in ("bed", "csv"):
        urls[f"human_liver_native.{fmt}"] = BASE + "download/export/?" + urlencode({"assembly": "GRCh38", "experiment": human_ids, "format": fmt})
    urls["mouse_liver_native.bed"] = BASE + "download/export/?" + urlencode({"assembly": "GRCm38", "experiment": "EXM0064", "format": "bed"})
    urls["mouse_tissue_native.bed"] = BASE + "download/export/?" + urlencode({"assembly": "GRCm38", "sample_class": "tissue", "target": "H3K27ac", "is_perturbed": "false", "format": "bed"})
    for req in json.loads((SRC / "comparator_download_requests.json").read_text()):
        urls[req["file"]] = req["url"]
        urls[req["file"].replace(".bed", ".count.json")] = req["url"].replace("/export/", "/summary/")
    save_json(ROOT / "source_manifest.json", [dict(
        file="sources/" + filename, url=url, retrieved_date="2026-11-02",
        bytes=(SRC / filename).stat().st_size,
        sha256=hashlib.sha256((SRC / filename).read_bytes()).hexdigest())
        for filename, url in sorted(urls.items())])
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
