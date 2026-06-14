import os
import re

import numpy as np

from libmeica.utils.selection import andb, dice

ENCODING_FSPACE_FILE = "encoding_fspace.txt"


class EncodingTable:
    """
    Column-oriented numeric table with attribute access.

    label_spec format:
        label: (numeric_width, default_value)

    dtype is inferred from default_value
    """

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------

    def __init__(self, *, Ne: int = 0, nc: int = 0):
        label_spec = {
            "N": (5, 0),
            "C": (5, 3),
            "G": (5, 3),
            "K": (5, 0.0),
            "R": (7, 0.0),
            "V": (7, 0.0),
            "score": (9, 0.0),
            "zu_KR": (9, 0.0),
            "zu_sr2_ss0": (11, 0.0),
            "zu05_R": (7, 0.0),
            "zu_er": (7, 0.0),
            "zu_ex": (7, 0.0),
            "ex": (7, 0),
            "sr2": (7, 0),
            "ss0": (7, 0),
            "sz": (7, 0),
            "tr2": (7, 0),
            "ts0": (7, 0),
            "tz": (7, 0),
            "lr2": (7, 0),
            "ls0": (7, 0),
            "lz": (7, 0),
            "emx": (7, 0.0),
            "emr": (7, 0.0),
            "er": (7, 0.0),
            "cT1": (9, 0.0),
            "cT2s": (9, 0.0),
            "zu_cT1": (11, 0.0),
            "zu_cT2s": (11, 0.0),
            "sr2_ss0": (7, 0),
            "K_R": (8, 0),
        }

        object.__setattr__(self, "label_spec", label_spec)
        object.__setattr__(self, "labels", list(label_spec.keys()))
        object.__setattr__(
            self, "_label_to_idx", {k: i for i, k in enumerate(self.labels)}
        )
        object.__setattr__(self, "nlabels", len(self.labels))
        object.__setattr__(self, "Ne", Ne)

        self._allocate(nc)

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------

    def _dtype_for(self, label):
        default = self.label_spec[label][1]
        if isinstance(default, bool):
            return int
        return type(default)

    def _allocate(self, nc):
        F = np.zeros((nc, self.nlabels), dtype=float)

        for label, idx in self._label_to_idx.items():
            _, default = self.label_spec[label]
            if default != 0:
                F[:, idx] = default

        object.__setattr__(self, "Fspace", F)

    def reset(self, new_nc):
        if self.Fspace.shape[0] != new_nc:
            self._allocate(new_nc)

    # ------------------------------------------------------------------
    # attribute access
    # ------------------------------------------------------------------

    def __getattr__(self, label):
        if label in ["Ne"]:
            return object.__getattribute__(self, label)

        if label not in self._label_to_idx:
            raise AttributeError(label)

        idx = self._label_to_idx[label]
        # dtype = self._dtype_for(label)
        return self.Fspace[:, idx]

    def __setattr__(self, label, value):
        # Allow normal assignment for internal attributes
        if label in {
            "label_spec",
            "labels",
            "_label_to_idx",
            "nlabels",
            "Fspace",
            "Ne",
        }:
            object.__setattr__(self, label, value)
            return

        # Only known column labels are assignable
        if label not in self._label_to_idx:
            raise AttributeError(label)

        idx = self._label_to_idx[label]
        arr = np.asarray(value)

        # Allocate ONLY if table is currently empty
        if self.Fspace.shape[0] == 0:
            self._allocate(arr.shape[0])
        else:
            # Otherwise, shapes must match
            if arr.shape[0] != self.Fspace.shape[0]:
                raise ValueError(
                    f"Length mismatch for column '{label}': "
                    f"{arr.shape[0]} != {self.Fspace.shape[0]}"
                )

        # Enforce integer semantics on assignment (not on access)
        dtype = self._dtype_for(label)
        if dtype is int:
            arr = arr.astype(int)

        # Assign into backing array
        self.Fspace[:, idx] = arr

    # ------------------------------------------------------------------
    # formatting
    # ------------------------------------------------------------------

    def fmt_for(self, label):
        width, default = self.label_spec[label]
        dtype = self._dtype_for(label)

        # leading space ensures visual separation and alignment
        if dtype is int:
            return f" %{width}d"
        else:
            precision = 1 if width < 6 else 2
            return f" %{width}.{precision}f"

    def column_widths(self):
        """
        Compute exact printed column widths based on real data
        and the actual format strings.
        """
        column_widths = {}

        for ii, label in enumerate(self.labels):
            fmt = self.fmt_for(label)

            # Determine the value that will print widest
            col = self.Fspace[:, ii]

            if col.size == 0:
                test_vals = [self.label_spec[label][1]]
            else:
                vmax = np.nanmax(col)
                vmin = np.nanmin(col)
                test_vals = [vmax, vmin]

            value_width = max(len(fmt % v) for v in test_vals)
            label_width = len(label)

            column_width = max(label_width, value_width) + 1

            if ii == 0:
                column_width = max(column_width - 3, label_width)

            column_widths[label] = column_width

        return column_widths

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    def save(self, fn=ENCODING_FSPACE_FILE, header=""):

        widths = self.column_widths()

        label_line = "".join(f"{label:>{widths[label]}}" for label in self.labels)

        sequence_line = f"Ne:{self.Ne}"
        full_header = "\n".join(filter(None, [header, sequence_line, label_line]))

        try:
            hostname = os.environ["HOSTNAME"]
        except:
            hostname = ""

        try:
            meica_startdir = os.environ["STARTDIR"]
        except:
            meica_startdir = os.path.abspath(os.curdir)

        footer = f"{meica_startdir} {hostname}"

        fmts = [self.fmt_for(l) for l in self.labels]

        np.savetxt(
            fn,
            self.Fspace,
            fmt=fmts,
            header=full_header,
            footer=footer,
        )

    def load(self, fn=ENCODING_FSPACE_FILE):
        with open(fn) as f:
            header_lines = [
                ln for il, ln in enumerate(f) if ln.startswith("#") and il < 5
            ]

        if not header_lines:
            raise ValueError("No header found in file")

        Ne = 0
        for ln in header_lines:
            if "Ne:" in ln:
                m = re.search(r"\bNe\s*:\s*(\d+)\b", ln)
                if m:
                    Ne = int(m.group(1))
        if Ne == 0:
            raise ValueError("No TE count in header")
        self.Ne = int(Ne)

        labels = header_lines[-1].lstrip("#").split()
        data = np.loadtxt(fn)

        if data.ndim == 1:
            data = data[:, None]

        self.reset(data.shape[0])

        for i, label in enumerate(labels):
            if label not in self._label_to_idx:
                continue
            self.Fspace[:, self._label_to_idx[label]] = data[:, i]

        self.Fspace = self.Fspace[np.argsort(self.Fspace[:, 0])]

    # ------------------------------------------------------------------
    # domain logic
    # ------------------------------------------------------------------

    def tdice(self, A, B):
        mask = self.C != 0
        return dice(A[mask], B[mask])

    def gtpenalty(self, G=None):
        if G is None:
            G = self.G
        acc_to_ej = self.tdice(G == 1, self.C == 3)
        acc_to_ign = self.tdice(G == 1, self.C == 2)
        acc_to_midk = self.tdice(G == 1, self.C == -1)
        acc_to_rej = dice(G == 1, self.C == 0)

        ej_to_acc = self.tdice(G == 2, self.C == 1)
        ign_to_acc = self.tdice(G == 3, self.C == 1)
        ign_to_acc_sel = andb([G == 3, self.C == 1]) == 2
        rej_to_acc = dice(G == 0, self.C == 1)
        midk_to_acc = self.tdice(G == -1, self.C == 1)

        return (
            10**8 * (acc_to_rej + rej_to_acc)
            + 10**5 * acc_to_midk
            + 10**4.5 * midk_to_acc
            + 10**4 * acc_to_ej
            + 10**3.5 * ej_to_acc
            + 10**3 * acc_to_ign
            + 10 * ign_to_acc * (self.V[ign_to_acc_sel].sum()) ** 2
        )
