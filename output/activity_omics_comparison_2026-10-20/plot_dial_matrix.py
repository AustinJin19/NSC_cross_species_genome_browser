#!/usr/bin/env python3
"""Dial correlogram of the existing, jointly matched MLS association profiles.

The filled fraction of each dial equals abs(Spearman rho). Positive dials fill
clockwise (red), negative dials counterclockwise (blue), from twelve o'clock.
The upper triangle reports the same coefficients numerically.
"""
import csv
import math
import os
import sys
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / 'cross_species_visual_tool/web_app').is_dir())
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'cross_species_visual_tool/web_app/.deps'))
os.environ.setdefault('MPLCONFIGDIR', str(HERE / 'mpl_cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Rectangle

OUT = ROOT / 'outputs/activity_omics_comparison_2026-10-20'
METHODS = ['tss2kb', 'tss5kb', 'body', 'tss200kb', 'tss500kb', 'rna', 'protein']
LABELS = ['TSS ±2 kb', 'TSS ±5 kb', 'Gene body + 2 kb', 'TSS ±200 kb',
          'TSS ±500 kb', 'Transcriptome', 'Proteome']
POSITIVE, NEGATIVE = '#BB3E3F', '#246A9C'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': .85})


def dial(ax, x, y, value, radius=.34):
    color = POSITIVE if value >= 0 else NEGATIVE
    ax.add_patch(Circle((x, y), radius, facecolor='#F7F8FA', edgecolor='#B7C1CB', linewidth=.8))
    # Cartesian y coordinates are used so that angle +90 points upward.
    end = 90 - value * 360
    if value > 0:
        start_angle, end_angle = end, 90
    else:
        start_angle, end_angle = 90, end
    if value != 0:
        ax.add_patch(Wedge((x, y), radius, start_angle, end_angle,
                           facecolor=color, edgecolor='none'))
        angle = math.radians(end)
        ax.plot([x, x+radius*math.cos(angle)], [y, y+radius*math.sin(angle)],
                color=color, linewidth=.95, solid_capstyle='round')
    ax.plot([x, x], [y+radius*.91, y+radius*1.09], color='#6B7783', linewidth=.7)


def main():
    with (OUT / 'joint_correlation_map.tsv').open() as h:
        records = list(csv.DictReader(h, delimiter='\t'))
    lookup = {(r['method_a'], r['method_b']): r for r in records}
    assert len(lookup) == 49
    values = np.array([[float(lookup[(a, b)]['spearman_rho']) for b in METHODS] for a in METHODS])
    assert np.allclose(values, values.T) and np.allclose(np.diag(values), 1)
    assert np.all(np.isfinite(values)) and np.max(np.abs(values)) <= 1
    assert {int(r['tested_genes']) for r in records} == {3092}

    fig = plt.figure(figsize=(10.1, 9.3), facecolor='white')
    ax = fig.add_axes([.24, .25, .66, .715])
    ax.set_aspect('equal')
    ax.set_xlim(-.5, 6.5); ax.set_ylim(-.5, 6.5)
    # Reverse row placement without inverting the drawing coordinate system.
    for i in range(7):
        y = 6-i
        for j in range(7):
            val = values[i, j]
            if i == j:
                ax.add_patch(Rectangle((j-.5, y-.5), 1, 1, facecolor='#F0F3F6', edgecolor='none'))
                ax.text(j, y, '1.00', ha='center', va='center', fontsize=11, color='#7B8792')
            elif i > j:
                dial(ax, j, y, val)
            else:
                ax.text(j, y, f'{val:.2f}', ha='center', va='center', fontsize=13,
                        color=POSITIVE if val >= 0 else NEGATIVE)
    ax.set_xticks(range(7), LABELS)
    ax.set_yticks(range(7), LABELS[::-1])
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    ax.tick_params(length=0, pad=9, labelsize=11)
    ax.set_xticks(np.arange(-.5, 7, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 7, 1), minor=True)
    ax.grid(which='minor', color='#DEE4EA', linewidth=.6)
    ax.tick_params(which='minor', length=0)
    for spine in ax.spines.values(): spine.set_color('#344453')

    # The legend shows both direction and size of the filled sector, with a
    # common scale. It does not encode significance or correlations of abundance.
    legend = fig.add_axes([.24, .025, .66, .105])
    legend.set_aspect('equal'); legend.set_xlim(-.5, 6.5); legend.set_ylim(-.55, .65)
    legend.axis('off')
    for x, val in [(0, -1), (1.25, -.5), (2.5, 0), (3.75, .5), (5, 1)]:
        dial(legend, x, .22, val, .19)
        legend.text(x, -.13, f'{val:+.1f}' if val else '0', ha='center', va='top', fontsize=10)
    fig.text(.57, .015, 'Spearman ρ between gene–MLS correlation profiles', ha='center', fontsize=10, color='#344453')
    for suffix in ('png', 'pdf', 'svg'):
        target = OUT / f'method_correlation_dial_matrix.{suffix}'
        fig.savefig(target, dpi=300, bbox_inches='tight', pad_inches=.15, facecolor='white')
        print(target)
    plt.close(fig)


if __name__ == '__main__':
    main()
