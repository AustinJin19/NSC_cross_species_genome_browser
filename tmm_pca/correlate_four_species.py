#!/usr/bin/env python3
"""Correlate species-median TMM log-CPM profiles across all supplied genes."""
import csv
import hashlib
import json
from itertools import combinations
from pathlib import Path
import run_pca_all_genes as base
from scipy.stats import spearmanr, rankdata

np, plt = base.np, base.plt
HERE = Path(__file__).resolve().parent
OUT = HERE / "four_species_all_genes/correlation"
ORDER = ["Cricetulus_griseus", "Chinchilla_lanigera", "Octodon_degus", "Rattus_norvegicus"]
NAMES = ["Chinese hamster", "Chinchilla", "Degu", "Rat"]


def write_table(name, fields, records):
    with (OUT/name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source = HERE/"logcpm_tmm_core.tsv"
    meta_path = HERE/"four_species_all_genes/sample_metadata.tsv"
    metadata = base.read_table(meta_path)
    with source.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        samples = next(reader)[1:]
        rows = list(reader)
    genes = [row[0] for row in rows]
    values = np.array([row[1:] for row in rows], dtype=float)
    assert len(genes) == len(set(genes)) and len(samples) == len(set(samples))
    assert values.shape == (6459, 75) and np.isfinite(values).all()
    sample_groups = {sp:[r["sample"] for r in metadata if r["species"]==sp] for sp in ORDER}
    assert [len(sample_groups[sp]) for sp in ORDER] == [3,2,2,3]
    # Aggregate on the supplied logarithmic scale; do not log the values again.
    profiles = np.array([np.median(values[:,[samples.index(s) for s in sample_groups[sp]]], axis=1) for sp in ORDER])
    pearson = np.corrcoef(profiles)
    spearman = spearmanr(profiles, axis=1).statistic
    np.testing.assert_allclose(spearman, np.corrcoef(rankdata(profiles, axis=1)), atol=1e-12)
    for matrix in (pearson, spearman):
        assert np.isfinite(matrix).all()
        np.testing.assert_allclose(matrix, matrix.T)
        np.testing.assert_allclose(np.diag(matrix), 1)
    for metric,matrix in [("pearson",pearson),("spearman",spearman)]:
        write_table(f"{metric}_correlation.tsv", ["species",*NAMES],
                    [{"species":name, **dict(zip(NAMES,matrix[i]))} for i,name in enumerate(NAMES)])
    pairs = [dict(species_1=NAMES[i], species_2=NAMES[j], genes=len(genes),
                  samples_1=len(sample_groups[ORDER[i]]), samples_2=len(sample_groups[ORDER[j]]),
                  pearson_r=pearson[i,j], spearman_rho=spearman[i,j]) for i,j in combinations(range(4),2)]
    write_table("pairwise_correlations.tsv", list(pairs[0]), pairs)
    write_table("species_median_logcpm.tsv", ["gene",*NAMES],
                [{"gene":gene, **dict(zip(NAMES,profiles[:,i]))} for i,gene in enumerate(genes)])
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":11, "pdf.fonttype":42})
    fig,axes = plt.subplots(1,2,figsize=(11.5,5.5))
    fig.subplots_adjust(left=.14,right=.90,bottom=.24,top=.79,wspace=.55)
    for ax,matrix,title in zip(axes,[pearson,spearman],["Pearson r","Spearman ρ"]):
        im=ax.imshow(matrix,cmap="RdBu_r",vmin=-1,vmax=1,interpolation="nearest")
        ax.set_xticks(range(4),NAMES,rotation=35,ha="right")
        ax.set_yticks(range(4),NAMES)
        ax.set_xticks(np.arange(-.5,4,1),minor=True)
        ax.set_yticks(np.arange(-.5,4,1),minor=True)
        ax.grid(which="minor",color="white",linewidth=2)
        ax.tick_params(which="both",length=0)
        ax.set_title(title,pad=12,fontsize=13)
        for i in range(4):
            for j in range(4):
                ax.text(j,i,f"{matrix[i,j]:.3f}",ha="center",va="center",fontsize=12,
                        color="white" if abs(matrix[i,j])>.78 else "#202020")
    cax=fig.add_axes([.93,.27,.018,.48])
    fig.colorbar(im,cax=cax,ticks=[-1,-.5,0,.5,1],label="Correlation")
    fig.suptitle("Between-species ATAC activity correlation",fontsize=15,y=.975)
    fig.text(.52,.90,f"All {len(genes):,} genes · species medians of supplied TMM log-CPM",ha="center",fontsize=11)
    for ext in ["png","pdf"]:
        fig.savefig(OUT/f"species_correlation_heatmap.{ext}",dpi=240,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    summary=dict(genes=len(genes),species=NAMES,sample_groups=sample_groups,
                 aggregation="median across samples on supplied TMM log-CPM scale",
                 correlated_dimension="genes",renormalized=False,batch_corrected=False,
                 within_species_centered=False,
                 source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                 metadata_sha256=hashlib.sha256(meta_path.read_bytes()).hexdigest(),pairwise=pairs)
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"METHODS.txt").write_text(
        "Input: ../../logcpm_tmm_core.tsv; samples from ../sample_metadata.tsv.\n"
        "All 6,459 input genes; Chinese hamster (3 samples), chinchilla (2), degu (2), rat (3).\n"
        "Each species is represented by the per-gene median across its available samples on the supplied TMM log-CPM scale.\n"
        "For two samples this is the arithmetic mean of their logged values, not the log of their mean raw CPM.\n"
        "No additional log transform, variance filtering, re-normalization, species centering or batch correction.\n"
        "Pearson correlation compares logged activity profiles; Spearman correlation compares their ranks with average ranks for ties.\n"
        "Every species pair uses the same 6,459 genes. These correlate activity across genes, not activity with MLS.\n"
        "Original TMM factors were estimated in the full sample set and are retained.\n"
        "No p-values are reported; the coefficients summarize profile similarity and do not establish absence of batch effects.\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
