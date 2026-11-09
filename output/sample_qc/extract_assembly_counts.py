"""Extract the QC samples' assembly counts from the supplied reference workbook.

Reads the workbook without changing it. Requires openpyxl.
Selects one reference per species by matching the local bundle manifest's FASTA
path; retained alternative assemblies are never added together.
"""
import csv
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
WORKBOOK = ROOT / "outputs/reference_inventory_20261016/NSC_genome_assemblies_and_annotations.xlsx"
SHEET = "Local file evidence"

with (ROOT / "bundle_sample_qc_summary.tsv").open() as handle:
    qc = list(csv.DictReader(handle, delimiter="\t"))
with (ROOT / "data/liver_atac_bundle/manifest.local.tsv").open() as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))

keys = {(r["species"], r["sample"], r["run"]) for r in qc}
matched = [r for r in manifest if (r["species"], r["sample"], r["run"]) in keys]
if len(matched) != len(qc):
    raise ValueError("Expected one manifest record for each QC sample")

workbook = openpyxl.load_workbook(WORKBOOK, read_only=True, data_only=True)
rows = list(workbook[SHEET].iter_rows(values_only=True))
header_index = next(i for i, row in enumerate(rows) if row[0] == "Scientific name")
headers = rows[header_index]
evidence = [
    (index + 1, dict(zip(headers, row)))
    for index, row in enumerate(rows)
    if index > header_index and row[0]
]
results = []
for species in dict.fromkeys(r["species"] for r in qc):
    fasta_paths = {r["fasta"] for r in matched if r["species"] == species}
    if len(fasta_paths) != 1:
        raise ValueError(f"More than one reference used for {species}: {fasta_paths}")
    fasta = next(iter(fasta_paths))
    selected = [
        (row_number, row) for row_number, row in evidence
        if row["Scientific name"] == species.replace("_", " ")
        and row["FASTA file"] == fasta
    ]
    if len(selected) != 1:
        raise ValueError(f"Expected one workbook row for {species}, {fasta}")
    row_number, row = selected[0]
    count = row["FASTA sequences"]
    if not isinstance(count, (int, float)) or count < 1 or int(count) != count:
        raise ValueError(f"Invalid sequence count for {species}: {count}")
    results.append({
        "species": species,
        "sequence_count": int(count),
        "fasta_file": fasta,
        "source_workbook": str(WORKBOOK.relative_to(ROOT)),
        "source_sheet": SHEET,
        "source_row": row_number,
    })
workbook.close()

target = ROOT / "assembly_sequence_counts.tsv"
with target.open("w") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(results)
print(f"Extracted {len(results)} assembly counts to {target}")
