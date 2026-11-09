#!/usr/bin/env python3
"""Sample correlations for all 32 species, using all supplied TMM log-CPM genes."""
import argparse
import csv
import hashlib
import json
import re
from itertools import combinations
from pathlib import Path
import run_pca_all_genes as base
from scipy.stats import spearmanr, rankdata
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import squareform
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Patch

np, plt = base.np, base.plt
HERE = Path(__file__).resolve().parent
OUT = HERE / "all_genes/sample_correlation"
BATCH_COLORS = {"batch1":"#2878b5", "batch2":"#f28e2b"}


def write_table(name, fields, rows):
    with (OUT/name).open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def natural_key(name):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)",name)]


def heatmap(matrix, metadata, species_order, names, metric, gene_count):
    n=len(metadata)
    fig=plt.figure(figsize=(25,25))
    ax=fig.add_axes([.28,.12,.63,.63])
    im=ax.imshow(matrix,cmap="RdBu_r",vmin=-1,vmax=1,interpolation="nearest")
    labels=[r["sample"] for r in metadata]
    ax.set_xticks(range(n),labels,rotation=90,ha="center",fontsize=8)
    ax.set_yticks(range(n),labels,fontsize=8)
    ax.tick_params(which="both",length=0,pad=3)
    boundaries=np.cumsum([sum(r["species"]==sp for r in metadata) for sp in species_order])
    for boundary in boundaries[:-1]:
        ax.axhline(boundary-.5,color="#555555",linewidth=.5,alpha=.85)
        ax.axvline(boundary-.5,color="#555555",linewidth=.5,alpha=.85)
    start=0
    for stop,sp in zip(boundaries,species_order):
        mid=(start+stop-1)/2
        # Common-name group labels are separate from the exact sample IDs.
        ax.text(-.102,mid,names[sp],transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=8.5)
        ax.plot([-.092,-.087,-.087,-.092],[start-.4,start-.4,stop-.6,stop-.6],
                transform=ax.get_yaxis_transform(),color="#666666",lw=.65,clip_on=False)
        ax.text(mid,1.026,names[sp],transform=ax.get_xaxis_transform(),rotation=90,
                ha="center",va="bottom",fontsize=8.5)
        start=stop
    for i,row in enumerate(metadata):
        ax.add_patch(Rectangle((i-.5,1.006),1,.009,transform=ax.get_xaxis_transform(),
                               color=BATCH_COLORS[row["seq_batch"]],linewidth=0,clip_on=False))
    cax=fig.add_axes([.935,.35,.013,.28])
    fig.colorbar(im,cax=cax,ticks=[-1,-.5,0,.5,1],label=metric)
    fig.text(.58,.925,f"{metric} · sample-to-sample ATAC activity correlation",ha="center",fontsize=20)
    fig.text(.58,.902,f"{n} samples · {len(species_order)} species · all {gene_count:,} genes",ha="center",fontsize=15)
    fig.legend(handles=[Patch(facecolor=c,label=b) for b,c in BATCH_COLORS.items()],
               loc="upper center",bbox_to_anchor=(.58,.888),ncol=2,frameon=False,
               title="Sequencing batch",fontsize=12,title_fontsize=12)
    return fig


def clustered_heatmap(matrix, metadata, tree, metric, gene_count):
    """Use one correlation-distance tree for both axes, with exact sample labels."""
    n = len(metadata)
    fig = plt.figure(figsize=(25, 25))
    ax = fig.add_axes([.21, .12, .65, .65])
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    labels = [row["sample"] for row in metadata]
    ax.set_xticks(range(n), labels, rotation=90, ha="center", fontsize=8)
    ax.set_yticks(range(n), labels, fontsize=8)
    ax.yaxis.tick_right()
    ax.tick_params(which="both", length=0, pad=3)

    top = fig.add_axes([.21, .79, .65, .12])
    left = fig.add_axes([.075, .12, .12, .65])
    for tree_ax, orientation in [(top, "top"), (left, "left")]:
        rendered = dendrogram(tree, ax=tree_ax, orientation=orientation,
                              no_labels=True, color_threshold=0,
                              above_threshold_color="#333333")
        assert rendered["leaves"] == leaves_list(tree).tolist()
        for collection in tree_ax.collections:
            collection.set_linewidth(.8)
        tree_ax.set_axis_off()
    left.invert_yaxis()

    batch_top = fig.add_axes([.21, .775, .65, .009])
    batch_left = fig.add_axes([.198, .12, .009, .65])
    for i, row in enumerate(metadata):
        color = BATCH_COLORS[row["seq_batch"]]
        batch_top.add_patch(Rectangle((i, 0), 1, 1, color=color, linewidth=0))
        batch_left.add_patch(Rectangle((0, i), 1, 1, color=color, linewidth=0))
    batch_top.set_xlim(0, n)
    batch_top.set_ylim(0, 1)
    batch_left.set_xlim(0, 1)
    batch_left.set_ylim(n, 0)
    batch_top.set_axis_off()
    batch_left.set_axis_off()

    cax = fig.add_axes([.945, .34, .013, .28])
    fig.colorbar(im, cax=cax, ticks=[-1, -.5, 0, .5, 1], label=metric)
    fig.text(.52, .975, f"{metric} · sample correlation clustered by similarity",
             ha="center", fontsize=20)
    fig.text(.52, .952, f"{n} samples · 32 species · all {gene_count:,} genes · average linkage, distance = 1 − correlation",
             ha="center", fontsize=13)
    fig.legend(handles=[Patch(facecolor=c, label=b) for b, c in BATCH_COLORS.items()],
               loc="upper center", bbox_to_anchor=(.52, .942), ncol=2,
               frameon=False, fontsize=12)
    return fig


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clustered", action="store_true",
                        help="Order samples by average-linkage clustering of 1 minus correlation.")
    args = parser.parse_args()
    if args.clustered:
        OUT = HERE / "all_genes/sample_correlation_clustered"
    OUT.mkdir(parents=True,exist_ok=True)
    source=HERE/"logcpm_tmm_core.tsv"
    meta_source=HERE/"pca_scores.tsv"
    original=base.read_table(meta_source)
    common_text=(HERE.parent/"plot_bundle_sample_qc.R").read_text()
    names=dict(re.findall(r'^  ([A-Za-z_]+) = "([^"]+)"',common_text,re.M))
    species_order=sorted({r["species"] for r in original},key=lambda sp:names[sp].lower())
    metadata=[dict(sample=r[""],species=r["species"],common_name=names[r["species"]],
                   seq_batch=r["seq_batch"]) for r in original]
    metadata.sort(key=lambda r:(species_order.index(r["species"]),natural_key(r["sample"])))
    samples=[r["sample"] for r in metadata]
    with source.open() as handle:
        reader=csv.reader(handle,delimiter="\t")
        columns=next(reader)[1:]
        rows=list(reader)
    genes=[r[0] for r in rows]
    values=np.array([r[1:] for r in rows],dtype=float)
    assert len(genes)==len(set(genes))==6459
    assert len(samples)==len(set(samples))==75 and len(species_order)==32
    assert set(samples)==set(columns)
    profiles=values[:,[columns.index(s) for s in samples]].T
    assert profiles.shape==(75,6459) and np.isfinite(profiles).all()
    factors={r[""]:r for r in base.read_table(HERE/"tmm_factors.tsv")}
    assert set(factors)==set(samples)
    for row in metadata:
        assert row["species"]==factors[row["sample"]]["species"]
        assert row["seq_batch"]==factors[row["sample"]]["seq_batch"]
    pearson=np.corrcoef(profiles)
    spearman=spearmanr(profiles,axis=1).statistic
    np.testing.assert_allclose(spearman,np.corrcoef(rankdata(profiles,axis=1)),atol=1e-12)
    matrices = {"pearson": pearson, "spearman": spearman}
    trees, orders = {}, {}
    for metric,matrix in matrices.items():
        assert np.isfinite(matrix).all()
        np.testing.assert_allclose(matrix,matrix.T)
        np.testing.assert_allclose(np.diag(matrix),1)
        order = np.arange(len(samples))
        if args.clustered:
            distance = np.clip(1 - (matrix + matrix.T) / 2, 0, 2)
            np.fill_diagonal(distance, 0)
            trees[metric] = linkage(squareform(distance), method="average", optimal_ordering=True)
            order = leaves_list(trees[metric])
            assert sorted(order.tolist()) == list(range(len(samples)))
            write_table(f"{metric}_linkage.tsv", ["left", "right", "distance", "sample_count"],
                        [dict(zip(["left", "right", "distance", "sample_count"], row))
                         for row in trees[metric]])
        orders[metric] = order
        ordered_samples = [samples[i] for i in order]
        ordered_matrix = matrix[np.ix_(order, order)]
        write_table(f"{metric}_sample_correlation.tsv",["sample",*ordered_samples],
                    [{"sample":s,**dict(zip(ordered_samples,ordered_matrix[i]))}
                     for i,s in enumerate(ordered_samples)])
        if args.clustered:
            write_table(f"{metric}_sample_order.tsv", ["order", *metadata[0]],
                        [dict(order=i+1, **metadata[j]) for i,j in enumerate(order)])
    if not args.clustered:
        write_table("sample_order.tsv",["order",*metadata[0]],
                    [dict(order=i+1,**row) for i,row in enumerate(metadata)])
    pairs=[dict(sample_1=samples[i],sample_2=samples[j],species_1=metadata[i]["common_name"],
                species_2=metadata[j]["common_name"],batch_1=metadata[i]["seq_batch"],
                batch_2=metadata[j]["seq_batch"],genes=len(genes),pearson_r=pearson[i,j],
                spearman_rho=spearman[i,j]) for i,j in combinations(range(len(samples)),2)]
    assert len(pairs)==2775
    write_table("pairwise_correlations.tsv",list(pairs[0]),pairs)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"pdf.fonttype":42})
    with PdfPages(OUT/"sample_correlations_all_species.pdf") as pdf:
        for key,matrix,title in [("pearson",pearson,"Pearson r"),("spearman",spearman,"Spearman ρ")]:
            if args.clustered:
                order = orders[key]
                fig = clustered_heatmap(matrix[np.ix_(order, order)],
                                        [metadata[i] for i in order], trees[key], title, len(genes))
            else:
                fig=heatmap(matrix,metadata,species_order,names,title,len(genes))
            fig.savefig(OUT/f"{key}_sample_correlation.png",dpi=240,bbox_inches="tight",facecolor="white")
            pdf.savefig(fig,bbox_inches="tight",facecolor="white")
            plt.close(fig)
    summary=dict(genes=len(genes),samples=len(samples),species=len(species_order),
                 species_order=[names[sp] for sp in species_order],
                 batch_counts={b:sum(r["seq_batch"]==b for r in metadata) for b in BATCH_COLORS},
                 sample_order=samples,unique_pairs=len(pairs),input_scale="supplied TMM log-CPM",
                 aggregation="none",batch_correction=False,
                 source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                 metadata_sha256=hashlib.sha256(meta_source.read_bytes()).hexdigest(),
                 pairwise_pearson_range=[min(r["pearson_r"] for r in pairs),max(r["pearson_r"] for r in pairs)],
                 pairwise_spearman_range=[min(r["spearman_rho"] for r in pairs),max(r["spearman_rho"] for r in pairs)])
    if args.clustered:
        summary.pop("species_order")
        summary.pop("sample_order")
        summary.update(clustering="average linkage; 1 - correlation distance; optimal leaf ordering",
                       axis_order="same row/column order within each heatmap; separate tree for each metric",
                       sample_orders={key: [samples[i] for i in order] for key,order in orders.items()})
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"METHODS.txt").write_text(
        "All 75 individual samples from all 32 species in ../../logcpm_tmm_core.tsv.\n"
        "All 6,459 genes are used for every pair; no new sample or gene exclusions.\n" +
        ("Both axes use average-linkage hierarchical clustering with distance = 1 - sample correlation and optimal leaf ordering.\n"
         "Pearson and Spearman are clustered independently; each heatmap uses the same tree/order for its two axes.\n"
         "Species and batch do not constrain clustering; they are metadata only.\n"
         if args.clustered else
         "Both axes group species alphabetically by common name, then order full sample IDs naturally within species.\n") +
        "Sample labels remain exact input IDs. Common names follow ../../../plot_bundle_sample_qc.R.\n"
        "Species and sequencing-batch metadata come from ../../pca_scores.tsv and are checked against ../../tmm_factors.tsv.\n"
        "Every cell compares two sample profiles of supplied TMM log-CPM values across genes.\n"
        "No aggregation across samples, extra logarithm, species centering, re-normalization or batch correction.\n"
        "Pearson correlates values; Spearman correlates ranks with average ranks for ties.\n" +
        ("Heatmaps share a fixed -1 to +1 color scale; dendrograms show sample similarity.\n"
         if args.clustered else
         "Heatmaps share a fixed -1 to +1 color scale and identical sample order; black lines delimit species.\n") +
        "Colored strips identify sequencing batch. Matrix TSVs and 2,775 unique pairs retain full precision.\n")
    print(json.dumps({k:summary[k] for k in ["genes","samples","species","batch_counts","unique_pairs","pairwise_pearson_range","pairwise_spearman_range"]},indent=2))


if __name__=="__main__":
    main()
