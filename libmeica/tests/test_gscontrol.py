"""Regression tests for libmeica.gscontrol.gscontrol_raw.

Pins the per-voxel Legendre detrend semantics of the spatial global
signal estimation. 4.0.2 removed a leftover ``[0]`` index that made
every voxel share the first masked voxel's fitted drift.
"""

import os
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np
from scipy.special import lpmv

import libmeica.gscontrol as gscontrol_mod
from libmeica.gscontrol import gscontrol_raw


class GscontrolTestBase(unittest.TestCase):
    dtrank = 4

    def _make_data(self, nx=4, ny=4, nz=4, nt=64, ne=1, seed=7):
        rng = np.random.default_rng(seed)
        t = np.linspace(-1, 1, nt)
        Lmix = np.array([lpmv(0, vv, t) for vv in range(self.dtrank)]).T
        nvox = nx * ny * nz
        offsets = rng.normal(100.0, 10.0, size=nvox)
        coefs = rng.normal(0.0, 5.0, size=(nvox, self.dtrank))
        ts = offsets[:, None] + coefs @ Lmix.T
        OCcatd = ts.reshape(nx, ny, nz, nt)
        catd = np.repeat(OCcatd[:, :, :, None, :], ne, axis=3)
        assets = SimpleNamespace(
            Ne=ne, head=None, aff=np.eye(4), catd=catd,
            tsshape=(nx, ny, nz, nt),
        )
        return OCcatd, assets, Lmix

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self._cwd = os.getcwd()
        os.chdir(self._tmp.name)
        self._niwrite = gscontrol_mod.niwrite
        gscontrol_mod.niwrite = lambda *a, **k: None

    def tearDown(self):
        os.chdir(self._cwd)
        gscontrol_mod.niwrite = self._niwrite
        self._tmp.cleanup()

    def _expected_t1mi(self, OCcatd, Lmix, mask):
        """Independent reference: per-voxel Legendre detrend, temporal min."""
        Gmu = OCcatd.mean(-1)
        dat = OCcatd[mask] - Gmu[mask][:, np.newaxis]
        sol = np.linalg.lstsq(Lmix, dat.T, rcond=None)
        detr = dat - np.dot(sol[0].T, Lmix.T)
        sphis = detr.min(1)
        sphis -= sphis.mean()
        return sphis


class TestGscontrolPerVoxelDetrend(GscontrolTestBase):
    def test_t1mi_matches_per_voxel_detrend(self):
        OCcatd, assets, Lmix = self._make_data()
        mask = np.ones(OCcatd.shape[:3], dtype=bool)
        gscontrol_raw(OCcatd=OCcatd, assets=assets, dtrank=self.dtrank)
        expected = self._expected_t1mi(OCcatd, Lmix, mask)
        np.testing.assert_allclose(assets.t1mi[mask], expected,
                                   rtol=1e-8, atol=1e-8)

    def test_t1mi_nondegenerate(self):
        OCcatd, assets, Lmix = self._make_data(seed=11)
        mask = np.ones(OCcatd.shape[:3], dtype=bool)
        gscontrol_raw(OCcatd=OCcatd, assets=assets, dtrank=self.dtrank)
        self.assertGreater(assets.t1mi[mask].std(), 0.0)

    def test_zero_mean_voxels_are_excluded(self):
        # gscontrol derives its mask from zero-mean voxels, so "masking"
        # means zeroing the time series of those voxels.
        OCcatd, assets, Lmix = self._make_data()
        OCcatd[2, :, :] = 0.0
        assets.catd = np.repeat(OCcatd[:, :, :, None, :], assets.Ne, axis=3)
        gscontrol_raw(OCcatd=OCcatd, assets=assets, dtrank=self.dtrank)
        np.testing.assert_array_equal(assets.t1mi[2, :, :], 0)
        self.assertTrue(np.all(np.isfinite(assets.t1mi[:2, :, :])))

    def test_output_shapes_and_finiteness(self):
        OCcatd, assets, Lmix = self._make_data(ne=2)
        out = gscontrol_raw(OCcatd=OCcatd, assets=assets, dtrank=self.dtrank)
        self.assertEqual(out.shape, (4, 4, 4, 2, 64))
        self.assertTrue(np.all(np.isfinite(out)))
        self.assertTrue(np.all(np.isfinite(assets.t1mi)))


if __name__ == "__main__":
    unittest.main()
