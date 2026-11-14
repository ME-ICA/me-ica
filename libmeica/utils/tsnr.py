"""tSNR (temporal SNR) map utility for ME-ICA.

Computes the temporal SNR of a 4D time series, per voxel: the ratio
of the temporal mean to the temporal standard deviation. 4.0.2 exports
the tSNR map of the denoised (hikts) time series as
``TED/tsnr_hikts.nii`` / ``<outprefix>_tsnr.nii``.
"""

import sys

import numpy as np

__all__ = ["tsnr_map", "tsnr_from_nifti"]


def tsnr_map(ts):
    """
    Compute the temporal SNR map of a 4D time series.

    Parameters
    ----------
    ts : array_like
        4D array of shape ``(nx, ny, nz, nt)``.

    Returns
    -------
    numpy.ndarray
        3D array of shape ``(nx, ny, nz)``: per-voxel
        ``mean(t) / std(t)`` over time (population standard
        deviation). Voxels with zero temporal standard deviation
        (constant time courses) are set to 0.
    """
    ts = np.asarray(ts, dtype=np.float64)
    if ts.ndim != 4:
        raise ValueError(
            "tsnr_map expects a 4D (nx, ny, nz, nt) array, got ndim=%d" % ts.ndim
        )
    mu = ts.mean(axis=3)
    sd = ts.std(axis=3)
    out = np.zeros(mu.shape, dtype=np.float64)
    nz = sd > 0
    out[nz] = mu[nz] / sd[nz]
    return out


def tsnr_from_nifti(in_path, out_path):
    """
    Compute the tSNR map of a 4D NIfTI time series and save it.

    Parameters
    ----------
    in_path : str
        Path to a 4D NIfTI file (e.g. ``TED/hik_ts_OC_T1c.nii``).
    out_path : str
        Path where the 3D tSNR map is written (float32, same affine
        and header as the input).
    """
    import nibabel as nib

    img = nib.load(in_path)
    ts = np.asarray(img.dataobj, dtype=np.float64)
    tmap = tsnr_map(ts)
    out = nib.Nifti1Image(tmap.astype(np.float32), img.affine, img.header)
    out.set_data_dtype(np.float32)
    nib.save(out, out_path)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("usage: tsnr.py <in_4d.nii> <out_3d.nii>\n")
        sys.exit(1)
    tsnr_from_nifti(sys.argv[1], sys.argv[2])
