"""Regression tests for the EJC group line in tedana.writect.

4.0.2: components that fall into none of the ACC/REJ/MID/IGN groups are
reported in a dedicated ``#EJC`` header line, so the group lines account
for every component in the table.
"""

import os
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np

from tedana import writect


def _run_writect(acc, rej, midk, ign, nc=10):
    cwd = os.getcwd()
    tmp = tempfile.TemporaryDirectory()
    os.chdir(tmp.name)
    try:
        comptable = np.zeros((nc, 5))
        for i in range(nc):
            comptable[i] = [i + 1, nc - i, 0.5, 10.0, 10.0]
        assets = SimpleNamespace(nt=100)
        writect(comptable, assets=assets, ctname="comp_table.txt",
                varexpl=42.0, classarr=[list(acc), list(rej),
                                        list(midk), list(ign)])
        with open("comp_table.txt") as f:
            return f.read()
    finally:
        os.chdir(cwd)
        tmp.cleanup()


def _group_line(text, tag):
    """Parse a ``#TAG idx,idx,...`` header line (empty when unassigned)."""
    for line in text.splitlines():
        if line.startswith("#" + tag + " "):
            head = line.split("\t")[0]
            parts = head.split(None, 1)
            if len(parts) < 2:
                return []
            idx = parts[1].strip()
            return [] if idx == "" else [int(x) for x in idx.split(",")]
    raise AssertionError("missing #%s line" % tag)


class TestCtabEjcLabel(unittest.TestCase):
    def test_ejc_lists_unassigned_components(self):
        text = _run_writect(acc=[1, 2, 3], rej=[4, 5], midk=[6], ign=[7])
        self.assertEqual(sorted(_group_line(text, "EJC")), [0, 8, 9])
        self.assertEqual(sorted(_group_line(text, "ACC")), [1, 2, 3])
        self.assertEqual(sorted(_group_line(text, "REJ")), [4, 5])
        self.assertEqual(sorted(_group_line(text, "MID")), [6])
        self.assertEqual(sorted(_group_line(text, "IGN")), [7])

    def test_ejc_empty_when_all_classified(self):
        text = _run_writect(acc=[1, 2, 3], rej=[4, 5], midk=[6],
                            ign=[0, 7, 8, 9])
        self.assertEqual(_group_line(text, "EJC"), [])

    def test_groups_partition_the_table(self):
        acc, rej, midk, ign = [1, 2, 3], [4, 5], [6], [7, 8]
        text = _run_writect(acc, rej, midk, ign, nc=10)
        all_groups = (sorted(_group_line(text, "ACC"))
                      + sorted(_group_line(text, "REJ"))
                      + sorted(_group_line(text, "MID"))
                      + sorted(_group_line(text, "IGN"))
                      + sorted(_group_line(text, "EJC")))
        self.assertEqual(sorted(all_groups), list(range(10)))


if __name__ == "__main__":
    unittest.main()
