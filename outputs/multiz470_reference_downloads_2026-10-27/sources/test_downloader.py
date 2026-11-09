"""Small local-file tests; never downloads a genome or uses the network."""
import csv
import gzip
import hashlib
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'download_multiz470_references.sh'


def run(script, *args, code=0):
    result = subprocess.run(['bash', str(script), *map(str, args)], capture_output=True, text=True)
    assert result.returncode == code, (result.returncode, result.stdout, result.stderr)
    return result


rows = list(csv.DictReader((ROOT / 'multiz470_reference_downloads.tsv').open(), delimiter='\t'))
embedded = run(SCRIPT, '--list').stdout.splitlines()
assert len(embedded) == len(rows) == 33
assert len({r['scientific_name'] for r in rows}) == 29
fields = ['msa_id', 'species_directory', 'filename', 'format', 'bytes', 'md5', 'download_url', 'status']
assert embedded == ['\t'.join(r[k] for k in fields) for r in rows]

with tempfile.TemporaryDirectory(prefix='multiz-download-test-') as temp:
    temp = Path(temp)
    dest = temp / 'output with spaces'
    result = run(SCRIPT, '--dry-run', '--out-dir', dest)
    assert len(result.stdout.splitlines()) == 32
    assert not dest.exists()
    result = run(SCRIPT, '--dry-run', '--only', 'hg38,mm10', '--out-dir', dest)
    assert len(result.stdout.splitlines()) == 2 and not dest.exists()
    run(SCRIPT, '--dry-run', '--only', 'not_a_real_assembly', code=1)
    run(SCRIPT, '--dry-run', '--only', 'HLratNor7', code=2)

    source = temp / 'fixture.fa.gz'
    source.write_bytes(gzip.compress(b'>chrTest\nACGTACGTACGTACGT\n', mtime=0))
    payload = source.read_bytes()
    md5 = hashlib.md5(payload).hexdigest()
    entry = '\t'.join(['testGenome', 'Test_species', 'fixture.fa.gz', 'fasta.gz', str(len(payload)), md5, source.as_uri(), 'available'])
    script = temp / 'test_download.sh'
    script.write_text(re.sub(r'(?<=\n)hg38\t.*?(?=\nMULTIZ_MANIFEST)', lambda _: entry, SCRIPT.read_text(), flags=re.S))
    run(script, '--out-dir', dest)
    target = dest / 'Test_species/testGenome/fixture.fa.gz'
    assert target.read_bytes() == payload
    assert 'Verified existing' in run(script, '--out-dir', dest).stdout

    resume = temp / 'resume'
    part = resume / 'Test_species/testGenome/fixture.fa.gz.part'
    part.parent.mkdir(parents=True)
    part.write_bytes(payload[:10])
    run(script, '--out-dir', resume)
    assert part.with_suffix('').read_bytes() == payload and not part.exists()

    target.write_bytes(b'corrupt')
    run(script, '--out-dir', dest, code=1)
    assert target.read_bytes() == b'corrupt'

    bad = temp / 'bad_md5.sh'
    bad.write_text(script.read_text().replace(md5, '0' * 32))
    bad_dest = temp / 'bad'
    run(bad, '--out-dir', bad_dest, code=1)
    assert not (bad_dest / 'Test_species/testGenome/fixture.fa.gz').exists()
    assert (bad_dest / 'Test_species/testGenome/fixture.fa.gz.part').exists()

print('PASS: manifest, dry-run, selection, unavailable IDs, download, resume, skip, corruption and MD5 failure')
