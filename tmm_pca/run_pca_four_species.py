#!/usr/bin/env python3
"""PCA of four requested species using all supplied TMM log-CPM genes.

Run: python3 tmm_pca/run_pca_four_species.py
Uses the existing 75-sample normalization; PCA is refitted to the 10 samples.
"""
import csv
import hashlib
import json
from pathlib import Path
import run_pca_all_genes as base

np, plt, Line2D = base.np, base.plt, base.Line2D
HERE = Path(__file__).resolve().parent
OUT = HERE / "four_species_all_genes"
SPECIES = {
    "Cricetulus_griseus": ("Chinese hamster", "#3C5488"),
    "Chinchilla_lanigera": ("Chinchilla", "#00A087"),
    "Octodon_degus": ("Degu", "#E64B35"),
    "Rattus_norvegicus": ("Rat", "#4DBBD5"),
}
BATCHES = {"batch1": ("#2878b5", "o"), "batch2": ("#f28e2b", "^")}


def write_table(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def labels(ax, scores, metadata):
    ax.figure.canvas.draw()
    points = ax.transData.transform(scores[:, :2])
    positions = points + [10, 12]
    names = [base.sample_label(row["sample"]) for row in metadata]
    widths = np.array([len(name)*7.5+6 for name in names])
    for _ in range(120):
        force = np.zeros_like(positions)
        for i in range(len(names)):
            for j in range(i):
                dx, dy = positions[i] - positions[j]
                if abs(dx) < (widths[i]+widths[j])/2+5 and abs(dy) < 17:
                    direction = 1 if dy >= 0 else -1
                    force[i, 1] += direction*1.5
                    force[j, 1] -= direction*1.5
        positions += force
        if np.max(np.abs(force)) < .01:
            break
    for i, name in enumerate(names):
        ax.annotate(name, scores[i, :2], xytext=ax.transData.inverted().transform(positions[i]),
                    fontsize=9, color="#202020", ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", color="#999999", lw=.5)
                    if np.linalg.norm(positions[i]-points[i]) > 23 else None)


def plot(scores, explained, metadata, title, stem, gene_count):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.7))
    fig.subplots_adjust(left=.075, right=.98, bottom=.21, top=.83, wspace=.25)
    for i, row in enumerate(metadata):
        sp, batch = row["species"], row["seq_batch"]
        axes[0].scatter(*scores[i, :2], c=SPECIES[sp][1], marker=BATCHES[batch][1],
                        s=85, edgecolor="#333333", linewidth=.65, zorder=3)
        axes[1].scatter(*scores[i, :2], c=BATCHES[batch][0], marker=BATCHES[batch][1],
                        s=85, edgecolor="#333333", linewidth=.65, zorder=3)
    for ax, color_label in zip(axes, ["Colored by species", "Colored by sequencing batch"]):
        base.style_axis(ax, explained)
        ax.margins(x=.24, y=.25)
        ax.set_title(color_label, fontsize=11, pad=10)
        labels(ax, scores, metadata)
    species_handles = [Line2D([], [], marker="s", linestyle="", color=color,
                       label=f"{name} (n = {sum(r['species']==sp for r in metadata)})", markersize=7)
                       for sp, (name,color) in SPECIES.items()]
    axes[0].legend(handles=species_handles, loc="upper center", bbox_to_anchor=(.5,-.17),
                   ncol=2, frameon=False, fontsize=9, columnspacing=1.5)
    batch_handles = [Line2D([], [], marker=marker, linestyle="", color=color,
                     label=f"{batch} (n = {sum(r['seq_batch']==batch for r in metadata)})", markersize=7)
                     for batch,(color,marker) in BATCHES.items()]
    axes[1].legend(handles=batch_handles, loc="upper center", bbox_to_anchor=(.5,-.17),
                   ncol=2, frameon=False, fontsize=9)
    fig.suptitle(f"{title}\nAll {gene_count:,} genes · {len(metadata)} samples", fontsize=13.5, y=.975)
    for ext in ["png", "pdf"]:
        fig.savefig(OUT / f"{stem}.{ext}", dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    source = HERE / "logcpm_tmm_core.tsv"
    with source.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        samples = next(reader)[1:]
        raw = list(reader)
    genes = [r[0] for r in raw]
    matrix = np.array([r[1:] for r in raw], dtype=float).T
    original = {r[""]:r for r in base.read_table(HERE / "pca_scores.tsv")}
    factors = {r[""]:r for r in base.read_table(HERE / "tmm_factors.tsv")}
    ids = [i for i,sample in enumerate(samples) if original[sample]["species"] in SPECIES]
    metadata = [dict(sample=samples[i], species=original[samples[i]]["species"],
                     common_name=SPECIES[original[samples[i]]["species"]][0],
                     seq_batch=original[samples[i]]["seq_batch"]) for i in ids]
    values = matrix[ids]
    assert len(ids) == 10 and values.shape == (10, len(genes))
    assert len(genes) == len(set(genes)) and np.isfinite(values).all()
    assert {row["species"] for row in metadata} == set(SPECIES)
    for row in metadata:
        assert factors[row["sample"]]["seq_batch"] == row["seq_batch"]
        assert factors[row["sample"]]["species"] == row["species"]
    scores, loadings, explained = base.pca(values)
    centered = values.copy()
    for sp in SPECIES:
        mask = np.array([r["species"] == sp for r in metadata])
        centered[mask] -= centered[mask].mean(axis=0)
        np.testing.assert_allclose(centered[mask].mean(axis=0), 0, atol=1e-12)
    within_scores, within_loadings, within_explained = base.pca(centered)
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10, "pdf.fonttype":42})
    plot(scores, explained, metadata, "PCA of Chinese hamster, chinchilla, degu and rat",
         "pca_four_species", len(genes))
    plot(within_scores, within_explained, metadata, "Within-species centered PCA of the four species",
         "pca_four_species_within_centered", len(genes))
    score_rows = []
    for i,row in enumerate(metadata):
        entry = dict(row)
        for pc in range(5):
            entry[f"PC{pc+1}"] = scores[i,pc]
            entry[f"within_PC{pc+1}"] = within_scores[i,pc]
        score_rows.append(entry)
    write_table("pca_scores.tsv", score_rows)
    write_table("sample_metadata.tsv", metadata)
    write_table("variance_explained.tsv", [dict(pca=name, pc=f"PC{i+1}", variance_fraction=v,
                   cumulative_variance_fraction=float(e[:i+1].sum()))
                for name,e in [("four_species",explained),("within_species",within_explained)]
                for i,v in enumerate(e)])
    write_table("gene_loadings.tsv", [dict(gene=g, **{f"PC{k+1}":loadings[i,k] for k in range(5)})
                                      for i,g in enumerate(genes)])
    sd = base.species_design(metadata)
    b = np.array([r["seq_batch"] == "batch2" for r in metadata], dtype=float)
    full = np.c_[sd, b]
    assoc = []
    for name, pcs in [("four_species",scores),("within_species",within_scores)]:
        for i in range(5):
            y = pcs[:,i]
            r0 = y-sd@np.linalg.lstsq(sd,y,rcond=None)[0]
            r1 = y-full@np.linalg.lstsq(full,y,rcond=None)[0]
            assoc.append(dict(pca=name, pc=f"PC{i+1}", R2_species=base.r_squared(y,sd),
                              R2_batch=base.r_squared(y,np.c_[np.ones(len(b)),b]),
                              R2_batch_given_species=max(0.,float(1-(r1@r1)/(r0@r0)))))
    write_table("pc_group_associations.tsv", assoc)
    summary = dict(genes=len(genes), samples=len(ids), species=4,
                   sample_counts={name:sum(r["species"]==sp for r in metadata) for sp,(name,_) in SPECIES.items()},
                   batch_counts={b:sum(r["seq_batch"]==b for r in metadata) for b in BATCHES},
                   variance_first5=explained[:5].tolist(), within_variance_first5=within_explained[:5].tolist(),
                   retained_original_TMM_normalization=True, refitted_PCA=True,
                   gene_scaling="center only, no unit-variance scaling", input_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    (OUT / "summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT / "METHODS.txt").write_text(
        "Subset the supplied TMM log-CPM matrix to Chinese hamster, chinchilla, degu and rat.\n"
        "All 6,459 supplied genes and all 10 available samples from these species are included.\n"
        "Chinese hamster: CCH1, CCH3, CCH5. Chinchilla: CH1, CH4. Degu: DG9, Degu5. Rat: Rat3, Rat4, Rat6.\n"
        "The original 75-sample TMM normalization is retained. TMM factors are not recalculated on this subset.\n"
        "PCA is recomputed within this subset; it is not a filtered display of the 75-sample PCA.\n"
        "Sample-by-gene matrix, centered per gene across the 10 samples, no variance filter or unit-variance scaling.\n"
        "Separate within-species diagnostic: subtract each species mean per gene, then recompute PCA using all genes.\n"
        "Within-species centering is not batch correction and should not be used for between-species MLS analyses.\n"
        "PC signs are oriented so the largest absolute loading is positive; signs are arbitrary.\n"
        "Seq_batch metadata come from ../pca_scores.tsv and are checked against ../tmm_factors.tsv.\n"
        "Group R-squared statistics are descriptive; batch-given-species is the fraction of residual variation reduced by adding batch.\n"
        "No samples are added from other files or excluded by sample QC. Original files remain unchanged.\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
