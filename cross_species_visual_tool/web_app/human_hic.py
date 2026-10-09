"""Bounded, local interval queries for the hg38 HepG2 4DN contact matrix."""
from pathlib import Path
from functools import lru_cache
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / '.deps'))
PATH = Path(__file__).resolve().parents[2] / '4DNFICSTCJQZ.hic'
RESOLUTIONS = (1000, 2000, 5000, 10000, 25000, 50000, 100000,
               250000, 500000, 1000000, 2500000, 5000000, 10000000)


@lru_cache(maxsize=24)
def query(chrom, lo, hi):
    """Input and returned intervals are 1-based inclusive; matrix bins are BED."""
    base = dict(label='HepG2 in situ Hi-C', accession='4DNFICSTCJQZ',
                assembly='hg38', available=False, cells=[])
    try:
        import hictkpy
        if not PATH.exists():
            return dict(base, error='Local Hi-C file unavailable')
        # Bound matrix allocation even for chromosome-wide browser windows.
        overview = hictkpy.File(str(PATH), 1000000)
        name = chrom.removeprefix('chr')
        size = overview.chromosomes().get(name)
        if size is None:
            return dict(base, error='Chromosome absent from Hi-C file')
        start, end = max(0, lo-1), min(size, hi)
        if end <= start:
            return dict(base, error='Outside Hi-C chromosome')
        resolution = next((r for r in RESOLUTIONS if (end-start)/r <= 160), RESOLUTIONS[-1])
        f = hictkpy.File(str(PATH), resolution)
        norm = 'KR' if 'KR' in f.avail_normalizations() else 'NONE'
        warning = None
        try:
            matrix = f.fetch(f'{name}:{start}-{end}', normalization=norm).to_numpy()
        except RuntimeError:
            if norm == 'NONE':
                raise
            norm = 'NONE'
            warning = 'KR unavailable for this region; showing raw counts'
            matrix = f.fetch(f'{name}:{start}-{end}').to_numpy()
        origin = start // resolution * resolution
        cells = []
        for i in range(matrix.shape[0]):
            for j in range(i, matrix.shape[1]):
                value = float(matrix[i, j])
                # Null denotes masked/invalid bins, never a measured zero.
                cells.append([i, j, value if np.isfinite(value) else None])
        finite = matrix[np.isfinite(matrix) & (matrix > 0)]
        cap = float(np.percentile(finite, 98)) if finite.size else 1.
        return dict(base, available=True, cells=cells, origin=origin,
                    resolution=resolution, normalization=norm, colorMax=cap,
                    bins=matrix.shape[0], warning=warning, chromSize=size)
    except (ImportError, OSError, RuntimeError, ValueError) as exc:
        return dict(base, error=f'Hi-C unavailable: {exc}')
