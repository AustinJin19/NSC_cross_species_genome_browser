#!/usr/bin/env python3
"""Prune the official TimeTree 6 mammal tree to the project's 32 PCA species plus human and rhesus macaque.

Run from any directory: python3 build_timetree.py
Requires Biopython and matplotlib. Does not estimate or modify divergence times.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import csv
import hashlib
import io
import json
import math
import os
import statistics

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(HERE / "sources" / "mpl_cache"))
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

SOURCE = HERE / "sources" / "TT6_Mammals.nwk"
METADATA = PROJECT / "tmm_pca/all_genes/sample_correlation/sample_order.tsv"
SOURCE_URL = "https://timetree.org/downloads/biodiversity/TT6_Mammals.nwk"
DOWNLOAD_PAGE = "https://timetree.org/downloads"
PREFIX = "timetree_34_species"

rows = list(csv.DictReader(METADATA.open(), delimiter="\t"))
common = {r["species"]: r["common_name"] for r in rows}
added_species = {"Homo_sapiens": "Human", "Macaca_mulatta": "Rhesus macaque"}
common.update(added_species)
sample_count = Counter(r["species"] for r in rows)
selected = set(common)
assert len(selected) == 34, f"Expected 34 species, found {len(selected)}"

source_tree = Phylo.read(SOURCE, "newick")
source_terminals = source_tree.get_terminals()
source_counts = Counter(t.name for t in source_terminals)
assert all(source_counts[n] == 1 for n in selected), "Missing or ambiguous tip match"
tips = {t.name: t for t in source_terminals if t.name in selected}
selected_root = source_tree.common_ancestor([tips[n] for n in sorted(selected)])

def prune(node):
    if node.is_terminal():
        return Clade(branch_length=node.branch_length, name=node.name) if node.name in selected else None
    children = [v for c in node.clades if (v := prune(c)) is not None]
    if not children:
        return None
    if len(children) == 1:
        child = children[0]
        child.branch_length = (child.branch_length or 0) + (node.branch_length or 0)
        return child
    return Clade(branch_length=node.branch_length, name=node.name, clades=children)

root = prune(selected_root)
root.branch_length = None  # Discard the stem above the selected species' MRCA.
tree = Tree(root=root, rooted=True)
tree.ladderize(reverse=True)
output_tips = tree.get_terminals()
assert {t.name for t in output_tips} == selected and len(output_tips) == 34
assert all((n.branch_length or 0) >= 0 for n in tree.find_clades())

# Check every species pair against the full, unmodified source tree.
distance_errors = []
pairwise = []
for a, b in combinations(sorted(selected), 2):
    original = source_tree.distance(tips[a], tips[b])
    current = tree.distance(a, b)
    distance_errors.append(abs(original - current))
    pairwise.append((a, b, current, current / 2))
assert max(distance_errors) < 1e-7, max(distance_errors)
depth = tree.depths()
tip_depths = [depth[t] for t in output_tips]
height = statistics.median(tip_depths)
assert max(tip_depths) - min(tip_depths) < 1e-7, "Source is not ultrametric"

Phylo.write(tree, HERE / f"{PREFIX}.nwk", "newick", format_branch_length="%.8f")
reread = Phylo.read(HERE / f"{PREFIX}.nwk", "newick")
assert len(reread.get_terminals()) == 34
assert max(abs(reread.distance(a, b) - d) for a, b, d, _ in pairwise) < 1e-7
(HERE / "species_for_timetree_upload.txt").write_text(
    "\n".join(n.replace("_", " ") for n in sorted(selected)) + "\n"
)

def order_for(name):
    if name in added_species:
        return "Primates"
    if name in {"Artibeus_jamaicensis", "Eonycteris_spelaea", "Eptesicus_fuscus", "Molossus_molossus"}:
        return "Chiroptera"
    if name == "Scalopus_aquaticus":
        return "Eulipotyphla"
    if name == "Sylvilagus_floridanus":
        return "Lagomorpha"
    return "Rodentia"

colors = {"Primates": "#E64B35", "Rodentia": "#008C75", "Lagomorpha": "#E18727", "Chiroptera": "#3C5488", "Eulipotyphla": "#98649C"}
with (HERE / "species_mapping.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["plot_order", "common_name", "scientific_name", "timetree_tip", "match", "order", "n_PCA_samples", "selection_source", "root_to_tip_Myr"])
    for i, t in enumerate(output_tips, 1):
        writer.writerow([i, common[t.name], t.name.replace("_", " "), t.name, "exact", order_for(t.name), sample_count[t.name], "User-added species" if t.name in added_species else "TMM/PCA metadata", f"{depth[t]:.8f}"])

with (HERE / "pairwise_distances.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["species_1", "species_2", "patristic_distance_Myr", "MRCA_age_Ma"])
    for a, b, d, age in pairwise:
        writer.writerow([a, b, f"{d:.8f}", f"{age:.8f}"])

nodes = list(tree.get_nonterminals(order="preorder"))
with (HERE / "node_ages.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["node_id", "age_Ma", "descendant_count", "descendant_species"])
    for i, n in enumerate(nodes, 1):
        writer.writerow([f"N{i:02d}", f"{height-depth[n]:.8f}", len(n.get_terminals()), ";".join(t.name for t in n.get_terminals())])

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(14.4, 14.0), facecolor="white")
ax = fig.add_axes([0.055, 0.105, 0.455, 0.79])
y = {t: i for i, t in enumerate(output_tips)}
def set_y(n):
    if n not in y:
        y[n] = (set_y(n.clades[0]) + set_y(n.clades[-1])) / 2
    return y[n]
set_y(root)
for n in tree.find_clades(order="preorder"):
    x = height - depth[n]
    if n.clades:
        ax.plot([x, x], [y[n.clades[0]], y[n.clades[-1]]], color="#293A43", lw=1.25, solid_capstyle="butt", zorder=3)
        for c in n.clades:
            ax.plot([x, height-depth[c]], [y[c], y[c]], color="#293A43", lw=1.25, solid_capstyle="butt", zorder=3)
    else:
        ax.plot(0, y[n], "o", markersize=4.5, color=colors[order_for(n.name)], clip_on=False, zorder=4)

ax.set_xlim(94, 0)
ax.set_ylim(33.9, -1.3)
ax.set_yticks([])
ax.set_xticks([80, 60, 40, 20, 0])
ax.set_xlabel("Divergence time (million years ago)", labelpad=13, fontsize=10.5)
ax.tick_params(axis="x", length=4, color="#5C6B73", labelsize=10)
ax.grid(axis="x", color="#E7ECEF", linewidth=0.7, zorder=0)
for side in ("top", "left", "right"):
    ax.spines[side].set_visible(False)
ax.spines["bottom"].set_color("#53646E")
transform = blended_transform_factory(ax.transAxes, ax.transData)
for t in output_tips:
    ax.text(1.025, y[t], common[t.name], transform=transform, ha="left", va="center", fontsize=10.7, color="#172D39", clip_on=False)
    ax.text(1.605, y[t], t.name.replace("_", " "), transform=transform, ha="left", va="center", fontsize=10.2, fontstyle="italic", color="#52616C", clip_on=False)
ax.text(1.025, -1.12, "COMMON NAME", transform=transform, fontsize=9, weight="bold", color="#677680", clip_on=False)
ax.text(1.605, -1.12, "SCIENTIFIC NAME", transform=transform, fontsize=9, weight="bold", color="#677680", clip_on=False)
fig.text(0.055, 0.957, "Evolutionary relationships of the 34 study species", fontsize=18, weight="bold", color="#152F40")
fig.text(0.055, 0.932, "TimeTree 6  |  32 PCA species + human and rhesus macaque  |  Time in millions of years", fontsize=10.5, color="#52616C")
handles = [Line2D([], [], marker="o", ls="", color=color, markersize=5, label=f"{order} ({sum(order_for(n)==order for n in selected)})") for order, color in colors.items()]
fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.965, 0.919), ncol=5, frameon=False, fontsize=9, handletextpad=0.4, columnspacing=1.3)
fig.text(0.055, 0.045, f"Source: official TimeTree 6 mammal tree, accessed 27 October 2026. Selected-species MRCA: {height:.2f} Ma.", fontsize=8.5, color="#52616C")
fig.text(0.055, 0.028, "All 34 species matched exactly; original topology and branch distances retained. Bowhead whale is not included.", fontsize=8.3, color="#52616C")
for extension in ("pdf", "svg", "png"):
    fig.savefig(HERE / f"{PREFIX}.{extension}", dpi=220, facecolor="white")
plt.close(fig)

manifest = {
    "source": SOURCE_URL, "download_page": DOWNLOAD_PAGE, "source_version": "TimeTree 6",
    "accessed": "2026-10-27", "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "source_tip_count": len(source_terminals), "selection_source": str(METADATA.relative_to(PROJECT)),
    "species_count": 34, "exact_species_matches": 34, "substitutions": [],
    "added_to_PCA_set": ["Homo_sapiens", "Macaca_mulatta"],
    "excluded_from_tree": ["Balaena_mysticetus"],
    "included_PCA_only_species": "Cricetomys_gambianus", "branch_length_unit": "million years",
    "MRCA_age_Ma": height, "root_to_tip_range_Myr": [min(tip_depths), max(tip_depths)],
    "pairwise_distance_checks": len(pairwise), "maximum_distance_error_Myr": max(distance_errors),
    "negative_branches": 0, "zero_length_branches": sum(n.branch_length == 0 for n in tree.find_clades() if n is not root),
    "method": "Extract selected MRCA, prune nonselected tips, and sum adjacent branches when collapsing unary nodes. Discard only the stem above the selected MRCA. No new dating, imputation or species substitution.",
    "citation": "Kumar S et al. (2022). TimeTree 5: An Expanded Resource for Species Divergence Times. Molecular Biology and Evolution 39(8):msac174. https://doi.org/10.1093/molbev/msac174",
    "citation_note": "The TimeTree 6 download page currently requests the 2022 TimeTree 5 resource citation.",
    "use": "Local personal research output. See https://timetree.org/downloads for source terms."
}
(HERE / "provenance_and_validation.json").write_text(json.dumps(manifest, indent=2)+"\n")
(HERE / "README.txt").write_text(
    "34-species TimeTree for the cross-species liver ATAC project\n"
    "=========================================================\n\n"
    "Species selection: the 32 species in the existing TMM/PCA sample metadata, plus human and rhesus macaque.\n"
    "Includes Gambian pouched rat, human and rhesus macaque; excludes bowhead whale.\n"
    "Every name is an exact match to a terminal in the official TimeTree 6 mammal tree.\n\n"
    "timetree_34_species.nwk: rooted, dated Newick tree with scientific names. Branch lengths are millions of years.\n"
    "timetree_34_species.pdf / .svg / .png: vector and raster figures.\n"
    "species_mapping.tsv: common names, scientific names, sample counts, and plotted tip order.\n"
    "species_for_timetree_upload.txt: scientific names, one per line, for the TimeTree website.\n"
    "node_ages.tsv: ages and descendant taxa for internal nodes.\n"
    "pairwise_distances.tsv: evolutionary path distances and species-pair MRCA ages.\n"
    "provenance_and_validation.json: source, processing method, citation, and numerical validation.\n"
    "build_timetree.py: reproducible pruning and figure-generation script.\n\n"
    "All 561 pairwise distances were compared against the source tree and preserved. The tree is ultrametric.\n"
    "This is a pruned published synthesis, not a newly inferred tree from the project's ATAC sequences.\n"
    "No confidence intervals are present in the source Newick. Published relationships may remain uncertain.\n\n"
    f"Source: {SOURCE_URL}\nDownloaded: 2026-10-27\n"
    "Resource citation requested by TimeTree: Kumar et al. (2022), doi:10.1093/molbev/msac174.\n"
    "Use: personal research/teaching; consult TimeTree terms before redistribution.\n"
)
print(json.dumps({k: manifest[k] for k in ["species_count", "exact_species_matches", "MRCA_age_Ma", "pairwise_distance_checks", "maximum_distance_error_Myr", "zero_length_branches"]}, indent=2))
