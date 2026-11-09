"""Read significant GO terms from the supplied workbook and draw a dot plot.

Extraction: bundled Python, this script --extract
Rendering: /usr/bin/python3, this script --plot
The source workbook and its enrichment statistics are never modified.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import textwrap

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parent / "Regulatory_500kb_MLS_and_GO_2026-10-19.xlsx"


def extract():
    import openpyxl
    workbook = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
    rows = workbook["GO terms"].iter_rows(values_only=True)
    for _ in range(5):
        next(rows)
    headers = next(rows)
    significant, counts = [], collections.Counter()
    for row_number, values in enumerate(rows, 7):
        row = dict(zip(headers, values))
        if not isinstance(row.get("GO BH FDR"), (float, int)):
            continue
        counts[row["Input gene list"]] += 1
        if row["GO BH FDR"] >= 0.05:
            continue
        k, n = row["Overlap (k)"], row["Mapped query (n)"]
        K, N = row["Background term (K)"], row["Background (N)"]
        assert 0 < k <= n and k <= K <= N
        assert math.isclose(row["Fold enrichment"], (k / n) / (K / N), rel_tol=1e-10)
        assert 0 < row["GO BH FDR"] <= 1
        row["Gene ratio"] = k / n
        row["Minus log10 GO FDR"] = -math.log10(row["GO BH FDR"])
        row["Workbook row"] = row_number
        significant.append(row)
    significant.sort(key=lambda x: (x["Input gene list"], x["GO BH FDR"], x["Raw p-value"], x["Workbook row"]))
    assert len({(r["Input gene list"], r["GO ID"]) for r in significant}) == len(significant)
    workbook.close()
    data = {
        "source": str(SOURCE), "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "sheet": "GO terms", "threshold": "GO BH FDR < 0.05",
        "all_terms_by_list": dict(counts),
        "significant_terms_by_list": dict(collections.Counter(r["Input gene list"] for r in significant)),
        "terms": significant,
    }
    (OUT / "significant_go_terms.json").write_text(json.dumps(data, indent=2) + "\n")
    with (OUT / "significant_go_terms.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(significant[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(significant)
    print(json.dumps({k: v for k, v in data.items() if k != "terms"}, indent=2))


def plot():
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/nsc_go_dotplot_mpl")
    try:
        import matplotlib
    except ModuleNotFoundError:
        project = Path(__file__).resolve().parents[3]
        sys.path.insert(0, str(project / "cross_species_visual_tool/web_app/.deps"))
        import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    from matplotlib.cm import ScalarMappable
    from matplotlib.ticker import MultipleLocator

    data = json.loads((OUT / "significant_go_terms.json").read_text())
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == data["source_sha256"]
    terms = [r for r in data["terms"] if r["Input gene list"] == "FDR positive"]
    assert len(terms) == 7 and all(r["Mapped query (n)"] == 69 for r in terms)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                         "axes.linewidth": 0.8, "axes.labelsize": 12,
                         "xtick.labelsize": 10, "ytick.labelsize": 11})
    fig = plt.figure(figsize=(13.2, 5.8), facecolor="white")
    gs = fig.add_gridspec(1, 1, left=0.40, right=0.81, bottom=0.14, top=0.79)
    cmap = LinearSegmentedColormap.from_list("significance", ["#4679B5", "#B5A2B6", "#BE324B"])
    norm = Normalize(vmin=-math.log10(0.05), vmax=max(r["Minus log10 GO FDR"] for r in terms))
    groups = [
        ("FDR positive", "Positive genes: gene FDR < 0.05", "71 genes; 69 mapped"),
    ]
    for index, (group, title, subtitle) in enumerate(groups):
        selected = [r for r in terms if r["Input gene list"] == group]
        ax = fig.add_subplot(gs[index])
        ax.set_facecolor("white")
        ax.scatter([r["Fold enrichment"] for r in selected], range(len(selected)),
                   s=[r["Overlap (k)"] * 11 for r in selected],
                   c=[r["Minus log10 GO FDR"] for r in selected], cmap=cmap, norm=norm,
                   edgecolors="#333333", linewidths=0.6, zorder=3)
        ax.set_yticks(range(len(selected)), [textwrap.fill(r["Biological process"], 43)
                                            for r in selected])
        ax.set_ylim(len(selected) - 0.45, -0.6)
        ax.set_xlim(0, 27)
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.set_axisbelow(True)
        ax.grid(axis="both", color="#E4E7EA", linewidth=0.7)
        ax.tick_params(axis="y", length=0, pad=12)
        ax.set_xlabel("Fold enrichment", labelpad=9)

    fig.text(0.035, 0.975, "GO Biological Process enrichment", fontsize=16,
             fontweight="bold", ha="left", va="top")
    fig.text(0.035, 0.916, "TSS ±500 kb regulatory activity–MLS associations", fontsize=11,
             color="#505860", ha="left", va="top")
    fig.text(0.035, 0.856, "71 positive genes (gene FDR < 0.05); 69 mapped · GO FDR < 0.05",
             fontsize=11, ha="left", va="top")
    cax = fig.add_axes([0.865, 0.43, 0.018, 0.35])
    bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    bar.set_label("−log₁₀(GO FDR)", labelpad=10, fontsize=11)
    bar.set_ticks([1.5, 2.5, 3.5, 4.5])
    legend_ax = fig.add_axes([0.845, 0.17, 0.13, 0.18])
    legend_ax.axis("off")
    handles = [legend_ax.scatter([], [], s=k * 11, c="#8E969E", edgecolors="#333333",
                                 linewidths=0.6, label=str(k)) for k in sorted({r["Overlap (k)"] for r in terms})]
    legend_ax.legend(handles=handles, title="Gene count", loc="upper left",
                     frameon=False, labelspacing=1.6, borderpad=0.5)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"significant_GO_dotplot_71genes.{extension}", dpi=300, facecolor="white")
    plt.close(fig)
    with (OUT / "significant_go_terms_71genes.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(terms[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(terms)
    (OUT / "METHODS_71genes.txt").write_text(
        "Source: " + str(SOURCE) + "\nSheet: GO terms; header row 6.\n"
        "All 7 terms with GO-term BH FDR < 0.05 from the 71 positive genes passing gene-level BH FDR < 0.05 (69 mapped) are plotted.\n"
        "Only the FDR positive gene list is shown; no enrichment was rerun.\n"
        "GO BH correction is separate per gene list, as recorded in the source.\n"
        "x = fold enrichment = (k/n)/(K/N); dot area is proportional to overlap k; color = -log10(GO FDR).\n"
        "Terms are sorted by GO FDR then raw p. Each term contains 7 query genes, so all dots have the same size.\n"
        "Source model: log10(species-median peak count) versus log10(MLS), adjusted for log10(body mass).\n"
        "TSS ±500 kb is a 1 Mb total window. Background: 18,032 mapped tested genes.\n"
        "Terms may overlap; plotting does not establish independent pathways or causality.\n"
    )
    print("Saved PNG, PDF and SVG:", OUT / "significant_GO_dotplot_71genes")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.extract:
        extract()
    if args.plot:
        plot()
