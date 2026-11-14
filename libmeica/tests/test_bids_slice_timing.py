"""Regression tests for libmeica.utils.bids.bids_slice_timing_cmd.

Pins the 4.0.2 behavior where the BIDS sidecar path occurs exactly once
in the generated extraction command. Pre-4.0.1 builds duplicated the
source directory in the embedded command, breaking ``--tpattern bids``.
"""

import json
import os
import subprocess
import tempfile
import unittest

from libmeica.utils.bids import bids_slice_timing_cmd


class TestBidsSliceTimingCmd(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self._cwd = os.getcwd()
        os.chdir(self._tmp.name)

    def tearDown(self):
        os.chdir(self._cwd)
        self._tmp.cleanup()

    def _sidecar(self, name="sub-01_task-rest_bold.json", timing=None):
        payload = {"RepetitionTime": 2.0}
        if timing is not None:
            payload["SliceTiming"] = timing
        path = os.path.join(self._tmp.name, name)
        with open(path, "w") as f:
            json.dump(payload, f)
        return path

    def test_path_occurs_exactly_once(self):
        path = self._sidecar(timing=[0.0, 0.5, 1.0, 1.5])
        st_1d, cmd = bids_slice_timing_cmd(path)
        self.assertEqual(cmd.count(f"open('{path}')"), 1)
        self.assertTrue(st_1d.endswith(".st.1D"))

    def test_command_writes_correct_st_1d(self):
        timing = [0.0, 0.25, 0.5, 0.75, 1.0]
        path = self._sidecar(timing=timing)
        st_1d, cmd = bids_slice_timing_cmd(path)
        subprocess.run(cmd, shell=True, check=True)
        with open(st_1d) as f:
            written = [float(x) for x in f.read().split()]
        self.assertEqual(written, timing)

    def test_missing_slice_timing_returns_none(self):
        path = self._sidecar(timing=None)
        self.assertIsNone(bids_slice_timing_cmd(path))


if __name__ == "__main__":
    unittest.main()
