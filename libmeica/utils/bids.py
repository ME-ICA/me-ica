"""BIDS sidecar helpers for ME-ICA."""

import json
import sys


def bids_slice_timing_cmd(json_path: str):
    """
    Build the shell command that extracts ``SliceTiming`` from a BIDS sidecar
    into a space-delimited ``.st.1D`` file for ``3dTshift -tpattern``.

    Parameters
    ----------
    json_path : str
        Full path to the BIDS sidecar (e.g. ``sub-01_task-rest_bold.json``).

    Returns
    -------
    (st_1d_path, cmd) : tuple
        Path of the ``.st.1D`` file the command writes, and the command
        string itself. The sidecar path occurs exactly once in the command.
    None
        When the sidecar carries no ``SliceTiming`` field.

    Notes
    -----
    Known limitation: the sidecar path is interpolated into a single-quoted
    ``python -c`` string, so paths containing single quotes will break.
    """
    with open(json_path) as ifh:
        sidecar = json.load(ifh)
    if "SliceTiming" not in sidecar:
        return None
    st_1d_path = json_path + ".st.1D"
    cmd = (
        f"{sys.executable} -c \"import json; open('{st_1d_path}', 'w').write("
        f"' '.join([str(tt) for tt in json.load(open('{json_path}'))['SliceTiming']]))\" "
    )
    return st_1d_path, cmd
