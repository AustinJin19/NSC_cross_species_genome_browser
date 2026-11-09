#!/usr/bin/env python3
"""Rerun global and within-species PCA using every supplied TMM log-CPM row.

Run: python3 tmm_pca/run_pca_all_genes.py
Requires numpy, scipy and matplotlib. Inputs are never modified.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEPS = ROOT / "cross_species_visual_tool/web_app/.deps"
if DEPS.exists():
    sys.path.insert(0, str(DEPS))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nsc_tmm_pca_mpl"))
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = HERE / "all_genes"


def read_table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(name, fields, records):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def pca(values):
    centered = values - values.mean(axis=0)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    # Deterministic orientation; the sign of a PC does not affect its meaning.
    signs = np.sign(vt[np.arange(len(vt)), np.abs(vt).argmax(axis=1)])
    signs[signs == 0] = 1
    scores = u * singular * signs
    loadings = vt.T * signs
    explained = singular ** 2 / np.sum(singular ** 2)
    np.testing.assert_allclose(scores @ loadings.T, centered, rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(explained.sum(), 1.0)
    return scores, loadings, explained


def r_squared(values, design):
    fitted = design @ np.linalg.lstsq(design, values, rcond=None)[0]
    return max(0.0, float(1 - np.sum((values-fitted)**2) / np.sum((values-values.mean())**2)))


def species_design(metadata):
    labels = sorted({r["species"] for r in metadata})
    return np.array([[r["species"] == sp for sp in labels] for r in metadata], dtype=float)


def associations(scores, explained, metadata, qc, tag):
    species = species_design(metadata)
    batches = sorted({r["seq_batch"] for r in metadata})
    batch = np.array([[r["seq_batch"] == b for b in batches] for r in metadata], dtype=float)
    full = np.c_[species, batch[:, 1:]]
    rows = []
    for i in range(5):
        y = scores[:, i]
        residual = y - species @ np.linalg.lstsq(species, y, rcond=None)[0]
        residual_full = y - full @ np.linalg.lstsq(full, y, rcond=None)[0]
        partial = max(0., float(1 - (residual_full @ residual_full) / (residual @ residual)))
        row = dict(pca=tag, pc=f"PC{i+1}", var_explained=explained[i],
                   R2_species=r_squared(y, species), R2_seq_batch=r_squared(y, batch),
                   R2_seq_batch_given_species=partial)
        for metric in ("tss_enrichment", "dup_pct", "final_read_pairs", "overall_align_pct"):
            values = np.array([qc[r["sample"]].get(metric, np.nan) for r in metadata])
            mask = np.isfinite(values)
            rho = float(spearmanr(y[mask], values[mask]).statistic)
            row[f"spearman_rho_{metric}"] = rho
            row[f"rho2_{metric}"] = rho*rho
            row[f"n_{metric}"] = int(mask.sum())
        rows.append(row)
    return rows


def style_axis(ax, explained, xpc=0, ypc=1):
    ax.set_facecolor("white")
    ax.grid(color="#e5e7eb", linewidth=0.55)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.8)
    ax.set_xlabel(f"PC{xpc+1} ({100*explained[xpc]:.2f}%)")
    ax.set_ylabel(f"PC{ypc+1} ({100*explained[ypc]:.2f}%)")
    ax.margins(x=.13, y=.15)


def sample_label(sample):
    return sample.removesuffix("Liver").replace("WtMice", "Mouse").replace("wtMice", "Mouse")


def add_labels(ax, x, y, names):
    # Offset nearby labels in display coordinates and retain leader lines.
    ax.figure.canvas.draw()
    points = ax.transData.transform(np.c_[x, y])
    positions = points + [5., 6.]
    widths = np.array([max(12., len(n)*4.2) for n in names])
    bounds = ax.get_window_extent()
    for iteration in range(160):
        shifts = np.zeros_like(positions)
        for i in range(len(names)):
            for j in range(i):
                dx, dy = positions[i]-positions[j]
                overlap_x = (widths[i]+widths[j])/2+3-abs(dx)
                overlap_y = 10.5-abs(dy)
                if overlap_x > 0 and overlap_y > 0:
                    force = min(overlap_y*.55, 3.0)
                    direction = 1 if dy >= 0 else -1
                    shifts[i, 1] += direction*force
                    shifts[j, 1] -= direction*force
                    if abs(dy) < 2:
                        shifts[i, 0] += .8 if dx >= 0 else -.8
                        shifts[j, 0] -= .8 if dx >= 0 else -.8
        positions += np.clip(shifts, -4, 4)
        positions[:, 0] = np.clip(positions[:, 0], bounds.x0+widths/2+3, bounds.x1-widths/2-3)
        positions[:, 1] = np.clip(positions[:, 1], bounds.y0+7, bounds.y1-7)
        if np.max(np.abs(shifts)) < .05:
            break
    for i, name in enumerate(names):
        xytext = ax.transData.inverted().transform(positions[i])
        ax.annotate(name, (x[i], y[i]), xytext=xytext, fontsize=5.6,
                    color="#303030", ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", color="#9ca3af", lw=.35)
                    if np.linalg.norm(positions[i]-points[i]) > 13 else None)


def save_figure(fig, name):
    for extension in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{extension}", dpi=240, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def plot_global(scores, explained, metadata, common):
    species = sorted({r["species"] for r in metadata}, key=lambda s: common.get(s, s))
    colors = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    color = dict(zip(species, colors))
    marker = {"batch1":"o", "batch2":"^"}
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 7.2))
    fig.subplots_adjust(left=.055, right=.775, bottom=.11, top=.89, wspace=.22)
    for ax, (a, b) in zip(axes, [(0,1), (2,3)]):
        for i, row in enumerate(metadata):
            ax.scatter(scores[i,a], scores[i,b], c=[color[row["species"]]],
                       marker=marker[row["seq_batch"]], s=44, edgecolor="#373737", lw=.5, zorder=3)
        style_axis(ax, explained, a, b)
        add_labels(ax, scores[:,a], scores[:,b], [sample_label(r["sample"]) for r in metadata])
    handles = [Line2D([], [], marker="s", linestyle="", color=color[sp],
                      label=common.get(sp, sp.replace("_", " ")), markersize=6) for sp in species]
    handles += [Line2D([], [], marker=m, linestyle="", color="#666666", label=batch, markersize=6)
                for batch, m in marker.items()]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(.785,.49),
               frameon=False, fontsize=7.3, labelspacing=.52)
    fig.suptitle(f"Global PCA · all {len(genes):,} consensus TSS genes · {len(metadata)} samples", fontsize=14, y=.96)
    save_figure(fig, "pca_global_all_genes")


def plot_within(scores, explained, metadata, qc):
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.6))
    fig.subplots_adjust(left=.07, right=.96, bottom=.12, top=.84, wspace=.27)
    colors = {"batch1":"#2878b5", "batch2":"#f28e2b"}
    markers = {"batch1":"o", "batch2":"^"}
    for batch in colors:
        ids = [i for i,r in enumerate(metadata) if r["seq_batch"]==batch]
        axes[0].scatter(scores[ids,0], scores[ids,1], label=batch, s=48,
                        c=colors[batch], marker=markers[batch], edgecolor="#373737", lw=.5, zorder=3)
    ids = [i for i,r in enumerate(metadata) if r["bridge_species"]=="yes"]
    axes[0].scatter(scores[ids,0], scores[ids,1], label="Species in both batches", s=110,
                    facecolors="none", edgecolors="#d62728", linewidths=1.0, zorder=4)
    axes[0].legend(loc="upper right", fontsize=8, frameon=True)
    axes[0].set_title("Sequencing batch", fontsize=11)
    tss = np.array([qc[r["sample"]]["tss_enrichment"] for r in metadata])
    finite = np.isfinite(tss)
    scatter = axes[1].scatter(scores[finite,0], scores[finite,1], c=tss[finite], cmap="viridis", s=48,
                              edgecolor="#373737", lw=.5, zorder=3)
    if (~finite).any():
        axes[1].scatter(scores[~finite,0], scores[~finite,1], c="#b8b8b8", s=48,
                        edgecolor="#373737", lw=.5, zorder=3, label="TSS enrichment unavailable")
        axes[1].legend(loc="upper right", fontsize=8, frameon=True)
    fig.colorbar(scatter, ax=axes[1], fraction=.04, pad=.025, label="TSS enrichment")
    axes[1].set_title("TSS enrichment", fontsize=11)
    for ax in axes:
        style_axis(ax, explained)
        add_labels(ax, scores[:,0], scores[:,1], [sample_label(r["sample"]) for r in metadata])
    fig.suptitle(f"Within-species centered PCA · all {len(genes):,} genes · {len(metadata)} samples", fontsize=14, y=.96)
    save_figure(fig, "pca_within_species_all_genes")


def plot_by_batch(scores, explained, metadata, gene_count):
    """Color the existing global PCA scores by sequencing batch."""
    colors = {"batch1": "#2878b5", "batch2": "#f28e2b"}
    markers = {"batch1": "o", "batch2": "^"}
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.6))
    fig.subplots_adjust(left=.07, right=.98, bottom=.12, top=.80, wspace=.23)
    for ax, (a, b) in zip(axes, [(0, 1), (2, 3)]):
        for batch in colors:
            ids = [i for i, row in enumerate(metadata) if row["seq_batch"] == batch]
            ax.scatter(scores[ids, a], scores[ids, b], c=colors[batch],
                       marker=markers[batch], s=52, edgecolor="#373737", lw=.5,
                       label=f"{batch} (n = {len(ids)})", zorder=3)
        style_axis(ax, explained, a, b)
        add_labels(ax, scores[:, a], scores[:, b],
                   [sample_label(row["sample"]) for row in metadata])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.53, .90),
               ncol=2, frameon=False, fontsize=10)
    fig.suptitle(f"All-gene PCA by sequencing batch · {gene_count:,} genes · {len(metadata)} samples",
                 fontsize=14, y=.97)
    save_figure(fig, "pca_global_by_batch_all_genes")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    matrix_path = HERE / "logcpm_tmm_core.tsv"
    with matrix_path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        samples = next(reader)[1:]
        records = list(reader)
    genes = [row[0] for row in records]
    values = np.array([row[1:] for row in records], dtype=float).T
    assert len(set(genes)) == len(genes) and len(set(samples)) == len(samples)
    assert values.shape == (len(samples), len(genes)) and np.isfinite(values).all()
    original = {r[""]: r for r in read_table(HERE / "pca_scores.tsv")}
    factors = {r[""]: r for r in read_table(HERE / "tmm_factors.tsv")}
    assert set(samples) == set(original) == set(factors)
    metadata = [{"sample":s, **{k:original[s][k] for k in ("species", "seq_batch", "bridge_species")}} for s in samples]
    for row in metadata:
        assert row["species"] == factors[row["sample"]]["species"]
        assert row["seq_batch"] == factors[row["sample"]]["seq_batch"]
    qc = {}
    for r in read_table(ROOT / "bundle_sample_qc_summary.tsv"):
        qc[r["sample"]] = {k:float(v) if v not in ("", "NA", "-") else np.nan
                            for k,v in r.items() if k not in ("sample", "species", "run")}
    assert all(s in qc for s in samples)
    common_text = (ROOT / "plot_bundle_sample_qc.R").read_text()
    common = dict(re.findall(r'^  ([A-Za-z_]+) = "([^"]+)"', common_text, re.M))
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":9, "pdf.fonttype":42})

    # Reproduce the old global PCA as a check on orientation and feature selection.
    top = np.argsort(np.var(values, axis=0, ddof=1))[-2000:]
    old_scores, _, _ = pca(values[:,top])
    reproduction = [abs(float(np.corrcoef(old_scores[:,i], [float(original[s][f"global_PC{i+1}"]) for s in samples])[0,1])) for i in range(5)]
    assert min(reproduction) > .99999

    global_scores, global_loadings, global_explained = pca(values)
    counts = {sp:sum(r["species"]==sp for r in metadata) for sp in {r["species"] for r in metadata}}
    keep = np.array([counts[r["species"]] >= 2 for r in metadata])
    within_meta = [r for r,k in zip(metadata,keep) if k]
    centered_species = values[keep].copy()
    for sp in counts:
        mask = np.array([r["species"]==sp for r in within_meta])
        if mask.any():
            centered_species[mask] -= centered_species[mask].mean(axis=0)
            np.testing.assert_allclose(centered_species[mask].mean(axis=0), 0, atol=1e-12)
    within_scores, within_loadings, within_explained = pca(centered_species)
    within_lookup = {r["sample"]:i for i,r in enumerate(within_meta)}
    output_scores = []
    for i,row in enumerate(metadata):
        entry = dict(row)
        for pc in range(5):
            entry[f"global_PC{pc+1}"] = global_scores[i,pc]
            entry[f"within_PC{pc+1}"] = within_scores[within_lookup[row["sample"]],pc] if row["sample"] in within_lookup else ""
        output_scores.append(entry)
    write_table("pca_scores.tsv", list(output_scores[0]), output_scores)
    for name, loadings in [("global",global_loadings),("within_species",within_loadings)]:
        write_table(f"{name}_gene_loadings.tsv", ["gene"]+[f"PC{i+1}" for i in range(5)],
                    [{"gene":gene, **{f"PC{i+1}":loadings[j,i] for i in range(5)}} for j,gene in enumerate(genes)])
    variance = [dict(pca=name, pc=f"PC{i+1}", variance_fraction=v, cumulative_variance_fraction=float(e[:i+1].sum()))
                for name,e in [("global",global_explained),("within_species",within_explained)] for i,v in enumerate(e)]
    write_table("variance_explained.tsv", list(variance[0]), variance)
    assoc = associations(global_scores, global_explained, metadata, qc, "global") + associations(within_scores, within_explained, within_meta, qc, "within_species")
    write_table("pc_covariate_assoc.tsv", list(assoc[0]), assoc)
    plot_global(global_scores, global_explained, metadata, common)
    plot_within(within_scores, within_explained, within_meta, qc)
    plot_by_batch(global_scores, global_explained, metadata, len(genes))
    summary = dict(genes=len(genes), global_samples=len(samples), species=len(counts),
                   within_species_samples=len(within_meta), within_species_species=sum(n>=2 for n in counts.values()),
                   global_variance_first5=global_explained[:5].tolist(),
                   within_species_variance_first5=within_explained[:5].tolist(),
                   constant_genes=int(np.sum(np.var(values,axis=0)==0)),
                   input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [matrix_path,HERE/"pca_scores.tsv",HERE/"tmm_factors.tsv",ROOT/"bundle_sample_qc_summary.tsv"]},
                   old_top2000_score_validation_abs_r=reproduction,
                   feature_selection="all supplied genes; no variance filter", gene_scaling="none; center only",
                   batch_correction=False, tmm_recomputed=False)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    (OUT / "METHODS.txt").write_text(
        f"Input: ../logcpm_tmm_core.tsv, all {len(genes):,} gene rows and {len(samples)} sample columns.\n"
        "Existing TMM-normalized logged values are used without further transformation, filtering or TMM recalculation.\n"
        "Global PCA: gene-centered sample-by-gene matrix; no unit-variance scaling; NumPy SVD.\n"
        f"Within-species PCA: {len(within_meta)} samples from {sum(n>=2 for n in counts.values())} species with at least two samples;\n"
        "subtract the species-specific mean of each gene before PCA. All supplied genes remain included.\n"
        "Within-species centering is a diagnostic operation, not a matrix for between-species MLS analysis.\n"
        "PC signs use the largest absolute loading as positive. Signs may differ from the earlier plots.\n"
        "Batch and species metadata come from ../pca_scores.tsv and are checked against ../tmm_factors.tsv.\n"
        "No samples are newly excluded by MLS, QC or batch. Bowhead, human and macaque are absent from the input matrix.\n"
        "QC values come from ../../bundle_sample_qc_summary.tsv, joined by exact sample ID.\n"
        "Samples without TSS enrichment remain in PCA and are plotted gray in the TSS color panel.\n"
        "Missing QC values are omitted only from the corresponding correlation, with its sample count exported.\n"
        "QC associations report signed Spearman correlations and their squares. Group associations use OLS R-squared.\n"
        "R2_seq_batch_given_species is the reduction in residual sum of squares after adding batch to a species model.\n"
        "FRiP and fraction-in-windows sample values were not supplied, so those associations and the FRiP color panel are omitted.\n"
        "Core library size divided by final read pairs is not substituted for an unavailable fraction-in-windows measurement.\n"
        "Original top-2,000 input files are preserved. Their first five global PCs were independently reproduced as validation.\n"
        "All components are listed in variance_explained.tsv; scores and gene loadings export the first five PCs.\n")
    print(json.dumps(summary, indent=2))
