"""Later signal downloads must survive regenerating the local bundle manifest."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import import_liver_bundle as bundle
import pyBigWig


class LaterDownloadTest(unittest.TestCase):
    def test_header_only_bigwig_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'species.fa.gz.fai').write_text('chr1\t1000\t0\t60\t61\n')
            with pyBigWig.open(str(root/'sample.cpm.bw'), 'w') as bw:
                bw.addHeader([('chr1',1000)])
            with self.assertRaisesRegex(ValueError, 'contains no signal records'):
                bundle.build('species', root, root)

    def test_new_signal_and_still_missing_sample(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            tracks = root/'data/atac/liver_atac_bundle/Ondatra_zibethicus'
            tracks.mkdir(parents=True)
            (tracks/'K11Liver.cpm.bw').write_bytes(b'signal')
            with patch.object(bundle, 'ROOT', root), patch.object(bundle, 'BUNDLE', root/'source'):
                found = bundle.manifest_track_path({'sample':'K11Liver','bigwig':'-'}, 'bigwig', tracks, {})
                self.assertEqual(found, 'data/atac/liver_atac_bundle/Ondatra_zibethicus/K11Liver.cpm.bw')
                self.assertEqual(bundle.manifest_track_path({'sample':'K12Liver','bigwig':'-'}, 'bigwig', tracks, {}), '-')

    def test_excluded_download_is_not_restored(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            tracks = root/'tracks'
            tracks.mkdir()
            (tracks/'NMR1Liver.cpm.bw').write_bytes(b'excluded')
            with patch.object(bundle, 'ROOT', root), patch.object(bundle, 'BUNDLE', root/'source'):
                result = bundle.manifest_track_path({'sample':'NMR1Liver','bigwig':'-'}, 'bigwig', tracks,
                                                    {'tracks/NMR1Liver.cpm.bw':None})
                self.assertEqual(result, '-')

    def test_ambiguous_signal_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for name in ('sample.bw', 'sample.cpm.bw'):
                (root/name).write_bytes(b'signal')
            with patch.object(bundle, 'ROOT', root):
                with self.assertRaises(ValueError):
                    bundle.manifest_track_path({'sample':'sample','bigwig':'-'}, 'bigwig', root, {})


if __name__ == '__main__':
    unittest.main()
