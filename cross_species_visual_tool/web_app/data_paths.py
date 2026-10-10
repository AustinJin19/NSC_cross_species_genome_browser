"""Locations of source workbooks; provenance continues to use original filenames."""
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[2]

def workbook_path(filename):
    folders = {'CrossSpecies_Count.xlsx': 'rna', 'rnaLongTPM.xlsx': 'rna',
               'protLongDF.xlsx': 'protein'}
    return PROJECT_ROOT / 'data' / 'omics' / folders[filename] / filename
