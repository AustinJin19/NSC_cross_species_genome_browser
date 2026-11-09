#!/usr/bin/env python3
"""Sample-level correlations across all genes, ordered by the four species."""
import csv
import hashlib
import json
import re
from itertools import combinations
from pathlib import Path
import run_pca_all_genes as base
from correlate_four_species import ORDER, NAMES
from scipy.stats import spearmanr, rankdata
from matplotlib.patches import Rectangle, Patch

np, plt = base.np, base.plt
HERE = Path(__file__).resolve().parent
OUT = HERE / "four_species_all_genes/sample_correlation"
BATCH_COLORS = {"batch1":"#2878b5", "batch2":"#f28e2b"}


def write_table(name, fields, rows):
    with (OUT/name).open("w", newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def natural_key(value):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)",value)]


def plot(pearson, spearman, metadata, gene_count):
    labels=[r["sample"] for r in metadata]
    n=len(labels)
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10, "pdf.fonttype":42})
    fig,axes=plt.subplots(1,2,figsize=(16.4,8.4))
    fig.subplots_adjust(left=.085,right=.915,bottom=.21,top=.76,wspace=.38)
    boundaries=np.cumsum([sum(r["species"]==sp for r in metadata) for sp in ORDER])
    for ax,matrix,title in zip(axes,[pearson,spearman],["Pearson r","Spearman ρ"]):
        im=ax.imshow(matrix,cmap="RdBu_r",vmin=-1,vmax=1,interpolation="nearest")
        ax.set_xticks(range(n),labels,rotation=55,ha="right",fontsize=9)
        ax.set_yticks(range(n),labels,fontsize=9)
        ax.set_xticks(np.arange(-.5,n,1),minor=True)
        ax.set_yticks(np.arange(-.5,n,1),minor=True)
        ax.grid(which="minor",color="white",linewidth=.65)
        ax.tick_params(which="both",length=0)
        for cut in boundaries[:-1]:
            ax.axhline(cut-.5,color="#343434",linewidth=1.3)
            ax.axvline(cut-.5,color="#343434",linewidth=1.3)
        start=0
        for stop,name in zip(boundaries,NAMES):
            ax.text((start+stop-1)/2,1.055,name.replace("Chinese hamster","Chinese\nhamster"),
                    transform=ax.get_xaxis_transform(),ha="center",va="bottom",fontsize=9)
            start=stop
        for i,row in enumerate(metadata):
            ax.add_patch(Rectangle((i-.5,1.008),1,.025,transform=ax.get_xaxis_transform(),
                                   color=BATCH_COLORS[row["seq_batch"]],linewidth=0,clip_on=False))
        ax.set_title(title,y=1.175,fontsize=14)
        for i in range(n):
            for j in range(n):
                ax.text(j,i,f"{matrix[i,j]:.3f}",ha="center",va="center",fontsize=8.4,
                        color="white" if abs(matrix[i,j])>.78 else "#202020")
    cax=fig.add_axes([.945,.28,.012,.40])
    fig.colorbar(im,cax=cax,ticks=[-1,-.5,0,.5,1],label="Correlation")
    fig.suptitle(f"Sample-to-sample ATAC correlation · all {gene_count:,} genes",fontsize=16,y=.99)
    handles=[Patch(facecolor=color,label=batch) for batch,color in BATCH_COLORS.items()]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.51,.945),ncol=2,
               frameon=False,fontsize=10,title="Sequencing batch")
    for ext in ("png","pdf"):
        fig.savefig(OUT/f"sample_correlation_heatmap.{ext}",dpi=260,bbox_inches="tight",facecolor="white")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    source=HERE/"logcpm_tmm_core.tsv"
    metadata=base.read_table(HERE/"four_species_all_genes/sample_metadata.tsv")
    metadata.sort(key=lambda r:(ORDER.index(r["species"]),natural_key(r["sample"])))
    samples=[r["sample"] for r in metadata]
    assert len(samples)==len(set(samples))==10
    with source.open() as handle:
        reader=csv.reader(handle,delimiter="\t")
        columns=next(reader)[1:]
        raw=list(reader)
    genes=[row[0] for row in raw]
    values=np.array([row[1:] for row in raw],dtype=float)
    profiles=values[:,[columns.index(s) for s in samples]].T
    assert profiles.shape==(10,6459) and np.isfinite(profiles).all()
    assert len(genes)==len(set(genes))
    pearson=np.corrcoef(profiles)
    spearman=spearmanr(profiles,axis=1).statistic
    np.testing.assert_allclose(spearman,np.corrcoef(rankdata(profiles,axis=1)),atol=1e-12)
    for metric,matrix in [("pearson",pearson),("spearman",spearman)]:
        np.testing.assert_allclose(matrix,matrix.T)
        np.testing.assert_allclose(np.diag(matrix),1)
        assert np.isfinite(matrix).all()
        write_table(f"{metric}_sample_correlation.tsv",["sample",*samples],
                    [{"sample":s,**dict(zip(samples,matrix[i]))} for i,s in enumerate(samples)])
    write_table("sample_order.tsv",["order",*metadata[0]],
                [dict(order=i+1,**row) for i,row in enumerate(metadata)])
    pairs=[dict(sample_1=samples[i],sample_2=samples[j],species_1=metadata[i]["common_name"],
                species_2=metadata[j]["common_name"],batch_1=metadata[i]["seq_batch"],
                batch_2=metadata[j]["seq_batch"],genes=len(genes),pearson_r=pearson[i,j],
                spearman_rho=spearman[i,j]) for i,j in combinations(range(len(samples)),2)]
    write_table("pairwise_correlations.tsv",list(pairs[0]),pairs)
    plot(pearson,spearman,metadata,len(genes))
    summary=dict(genes=len(genes),samples=len(samples),species_order=NAMES,sample_order=samples,
                 input_scale="supplied TMM log-CPM",aggregation="none",batch_correction=False,
                 input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                 pairwise_pearson_range=[min(r["pearson_r"] for r in pairs),max(r["pearson_r"] for r in pairs)],
                 pairwise_spearman_range=[min(r["spearman_rho"] for r in pairs),max(r["spearman_rho"] for r in pairs)])
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"METHODS.txt").write_text(
        "Compare the 10 individual samples from Chinese hamster, chinchilla, degu and rat.\n"
        "Input: ../../logcpm_tmm_core.tsv; all 6,459 supplied genes for every pair.\n"
        "Rows and columns have the identical species order, with natural sample-name ordering within species.\n"
        "Every cell compares two sample profiles of supplied TMM log-CPM values across genes.\n"
        "No averaging across samples, additional log transform, centering, re-normalization or batch correction.\n"
        "Pearson compares values; Spearman compares ranks, with average ranks for ties.\n"
        "Exact sample IDs label both axes. Group separators and headings mark species; top strips mark sequencing batch.\n"
        "The matrices and 45 unique sample pairs are exported at full precision.\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
