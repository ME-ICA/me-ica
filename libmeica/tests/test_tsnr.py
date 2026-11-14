"""Regression tests for libmeica.utils.tsnr (4.0.2 tSNR map export).

The 4.0.2 release computes the temporal-SNR map of the denoised
(hikts) time series and exports it as ``TED/tsnr_hikts.nii`` /
``<outprefix>_tsnr.nii``.
"""

import os
import tempfile
import unittest

import numpy as np

from libmeica.utils.tsnr import tsnr_from_nifti, tsnr_map


class TestTsnrMap(unittest.TestCase):
    def test_known_values(self):
        ts = np.zeros((1, 1, 2, 4))
        ts[0, 0, 0] = 10.0  # constant time course -> zero std -> tSNR 0
        ts[0, 0, 1] = [10.0, 20.0, 30.0, 40.0]
        out = tsnr_map(ts)
        self.assertEqual(out.shape, (1, 1, 2))
        self.assertEqual(out[0, 0, 0], 0.0)
        np.testing.assert_allclose(
            out[0, 0, 1], 25.0 / np.std([10.0, 20.0, 30.0, 40.0])
        )

    def test_matches_manual_mean_over_std(self):
        rng = np.random.default_rng(3)
        ts = rng.normal(100.0, 10.0, size=(3, 4, 5, 32))
        out = tsnr_map(ts)
        ref = ts.mean(axis=3) / ts.std(axis=3)
        np.testing.assert_allclose(out, ref, rtol=1e-12)

    def test_rejects_non_4d(self):
        with self.assertRaises(ValueError):
            tsnr_map(np.zeros((4, 5, 6)))

    def test_zero_volume_is_zero(self):
        out = tsnr_map(np.zeros((2, 2, 2, 8)))
        np.testing.assert_array_equal(out, 0.0)


class TestTsnrFromNifti(unittest.TestCase):
    def test_roundtrip(self):
        import nibabel as nib

        with tempfile.TemporaryDirectory() as tmp:
            rng = np.random.default_rng(11)
            ts = rng.normal(50.0, 5.0, size=(3, 3, 3, 16)).astype(np.float32)
            aff = np.diag([2.0, 2.0, 2.0, 1.0])
            in_path = os.path.join(tmp, "in.nii.gz")
            out_path = os.path.join(tmp, "tsnr.nii")
            nib.save(nib.Nifti1Image(ts, aff), in_path)
            tsnr_from_nifti(in_path, out_path)
            img = nib.load(out_path)
            self.assertEqual(img.shape, (3, 3, 3))
            np.testing.assert_array_equal(img.affine, aff)
            ref = tsnr_map(np.asarray(ts, dtype=np.float64))
            np.testing.assert_allclose(
                np.asarray(img.dataobj), ref, rtol=1e-5, atol=1e-5
            )


if __name__ == "__main__":
    unittest.main()
