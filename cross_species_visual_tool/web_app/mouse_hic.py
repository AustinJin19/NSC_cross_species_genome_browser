"""Mouse liver ZT10 contact matrix on the mm10 assembly."""
from functools import lru_cache
from pathlib import Path
from human_hic import query_matrix

PATH = Path(__file__).resolve().parents[2] / 'data/interactions/mouse_hic/mouse_liver_ZT10_merged.mcool'


@lru_cache(maxsize=24)
def query(chrom, lo, hi):
    base = dict(label='Mouse liver ZT10 Hi-C', assembly='mm10',
                available=False, cells=[])
    return query_matrix(PATH, chrom, lo, hi, base)
