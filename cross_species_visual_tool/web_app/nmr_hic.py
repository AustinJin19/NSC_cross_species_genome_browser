"""NMR hypothalamus contact matrix on the mHetGlaV3 assembly."""
from functools import lru_cache
from pathlib import Path
from human_hic import query_matrix

PATH = Path(__file__).resolve().parents[2] / 'data/interactions/nmr_microc/NMR_hypothalamus_merged.mcool'


@lru_cache(maxsize=24)
def query(chrom, lo, hi):
    base = dict(label='NMR hypothalamus Hi-C', assembly='mHetGlaV3',
                available=False, cells=[])
    return query_matrix(PATH, chrom, lo, hi, base)
