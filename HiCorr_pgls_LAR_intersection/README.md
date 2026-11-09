# HiCorr interactions overlapping pgls_LARs

Both input coordinate assemblies are hg19; the HiCorr assembly was confirmed by the user.

Inputs:
- `../Merged.raw_HiCorr_DeepLoop.1.5.bedpe`
- `../pgls_LARs_hg19.bed`

Output: `HiCorr_pgls_LARs_hg19.bedpe`

Retains a HiCorr interaction when **either anchor** overlaps a pgls_LAR by at least
one base. Original HiCorr columns 1–8, both full anchor intervals, row order, IDs,
and scores are preserved. No reciprocal overlap threshold or proximity window
is applied, and anchors are not trimmed to the intersection.

Added columns:

| Columns | Meaning |
|---|---|
| 9 | `anchor1`, `anchor2`, or `both` indicates which ends overlap LARs |
| 10–13 | Anchor 1: LAR count, IDs, hg19 intervals, and individual overlap lengths in bp |
| 14–17 | Anchor 2: LAR count, IDs, hg19 intervals, and individual overlap lengths in bp |

Multiple annotations are separated by semicolons and aligned positionally.
Unannotated anchors have count 0 and `.` in the other annotation fields.
Each input interaction is written once; multiple overlapping LARs do not duplicate
that interaction. Duplicate interactions already present in the input are preserved.
The header begins with `#` so BED-compatible tools can ignore it.

Overlap calculations use BED/BEDPE **0-based, half-open** coordinates. Touching
boundaries alone do not count. Added hg19 interval strings use that same convention.
The LAR ID is copied verbatim from BED column 4; some IDs contain coordinates
from a different assembly or convention. Overlap is based only on BED columns
1–3, never on coordinates embedded in the ID.

The original LAR BED score column is constant 1 and is not copied into the output.
Counts are in `summary.json`. Reproduce with `python3 intersect_hicorr.py`.
