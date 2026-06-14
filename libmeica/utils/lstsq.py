import multiprocessing as mp
import os
import platform
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from contextlib import contextmanager
from multiprocessing import shared_memory

import numpy as np
from scipy.linalg import solve_triangular

from .memory import trace_and_profile_inputs
from .volume import fmask, unmask


@contextmanager
def _temporary_env_vars(updates):
    """
    Temporarily set environment variables and restore them after.
    """
    original = {}
    try:
        for key, value in updates.items():
            original[key] = os.environ.get(key)
            os.environ[key] = str(value)
        yield
    finally:
        for key, value in original.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _worker_voxel_block(args):
    (
        start,
        end,
        shm_names,
        shapes,
        nb,
    ) = args

    shm_mdata = shared_memory.SharedMemory(name=shm_names["mdata"])
    shm_Q = shared_memory.SharedMemory(name=shm_names["Q"])
    shm_R = shared_memory.SharedMemory(name=shm_names["R"])
    shm_out = shared_memory.SharedMemory(name=shm_names["out"])

    mdata = np.ndarray(shapes["mdata"], dtype=np.float32, buffer=shm_mdata.buf)
    Q = np.ndarray(shapes["Q"], dtype=np.float32, buffer=shm_Q.buf)
    R = np.ndarray(shapes["R"], dtype=np.float32, buffer=shm_R.buf)
    out = np.ndarray(shapes["out"], dtype=np.float32, buffer=shm_out.buf)

    mchunk = mdata[:, start:end]

    # print("Worker starting", start, end)
    Qt_y = Q.T @ mchunk
    beta = solve_triangular(R, Qt_y, lower=False)

    if nb != 0:
        beta = beta[:-nb, :]

    out[start:end, :] = beta.T

    shm_mdata.close()
    shm_Q.close()
    shm_R.close()
    shm_out.close()


def _get_executor(n_jobs):
    return ThreadPoolExecutor(max_workers=n_jobs)
    # if sys.gettrace() is not None:
    #     # Debug mode
    #     return ThreadPoolExecutor(max_workers=n_jobs)
    # else:
    #     ctx = mp.get_context("spawn")
    #     return ProcessPoolExecutor(max_workers=n_jobs, mp_context=ctx)


def _get_n_jobs():
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        val = os.environ.get(var)
        if val:
            try:
                n_jobs = int(val)
                break
            except ValueError:
                pass
    else:
        n_jobs = os.cpu_count()
    return n_jobs


@trace_and_profile_inputs
def _get_coeffs_para(
    data,
    mask,
    X,
    *,
    add_const=False,
    betas5d=None,
    Ne5d=None,
    n_jobs=None,
):

    print("Running para get_coeffs")
    #
    # import pudb
    #
    # pudb.set_trace()

    if n_jobs is None:
        n_jobs = _get_n_jobs()

    # ---- Mask directly into shared memory buffer ----
    masked = fmask(data, mask).transpose().astype(np.float32)
    nt, Nm = masked.shape

    # ---- Design matrix ----
    X = np.atleast_2d(X)
    if X.shape[0] == 1:
        X = X.T

    if add_const is True:
        basemat = np.ones((X.shape[0], 1), dtype=X.dtype)
        X = np.hstack([X, basemat])
        nb = 1
    else:
        nb = 0

    Q, R = np.linalg.qr(X, mode="reduced")
    Q = np.ascontiguousarray(Q.astype(np.float32))
    R = np.ascontiguousarray(R.astype(np.float32))

    nc = Q.shape[1] - nb

    # ---- Shared memory allocations ----
    shm_mdata = shared_memory.SharedMemory(create=True, size=masked.nbytes)
    shm_Q = shared_memory.SharedMemory(create=True, size=Q.nbytes)
    shm_R = shared_memory.SharedMemory(create=True, size=R.nbytes)
    shm_out = shared_memory.SharedMemory(create=True, size=Nm * nc * 4)

    mdata_sh = np.ndarray((nt, Nm), dtype=np.float32, buffer=shm_mdata.buf)
    Q_sh = np.ndarray(Q.shape, dtype=np.float32, buffer=shm_Q.buf)
    R_sh = np.ndarray(R.shape, dtype=np.float32, buffer=shm_R.buf)
    out_sh = np.ndarray((Nm, nc), dtype=np.float32, buffer=shm_out.buf)

    mdata_sh[:] = masked
    Q_sh[:] = Q
    R_sh[:] = R

    # ---- Chunking ----
    # Use 2x cores for better load balance
    n_chunks = n_jobs * 2
    chunk_size = (Nm + n_chunks - 1) // n_chunks

    tasks = []
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, Nm)
        if start >= Nm:
            break

        tasks.append(
            (
                start,
                end,
                {
                    "mdata": shm_mdata.name,
                    "Q": shm_Q.name,
                    "R": shm_R.name,
                    "out": shm_out.name,
                },
                {
                    "mdata": (nt, Nm),
                    "Q": Q.shape,
                    "R": R.shape,
                    "out": (Nm, nc),
                },
                nb,
            )
        )

    # _ensure_linux_fork()

    # ---- Parallel ----
    with _temporary_env_vars(
        {
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
    ):
        with _get_executor(n_jobs) as ex:
            list(ex.map(_worker_voxel_block, tasks))

    tmpbetas = out_sh.copy()

    # ---- Cleanup ----
    shm_mdata.close()
    shm_mdata.unlink()
    shm_Q.close()
    shm_Q.unlink()
    shm_R.close()
    shm_R.unlink()
    shm_out.close()
    shm_out.unlink()

    # ---- Output ----
    if betas5d is not None:
        assert Ne5d is not None
        betas5d[:, :, :, Ne5d, :] = unmask(tmpbetas, mask)
    else:
        return unmask(tmpbetas, mask)


@trace_and_profile_inputs
def _get_coeffs_uni(data, mask, X, *, add_const=False, betas5d=None, Ne5d=None):
    """
    get_coeffs(data,X)

    Input:

    data has shape (nx,ny,nz,nt)
    mask has shape (nx,ny,nz)
    X    has shape (nt,nc)

    Output:

    out  has shape (nx,ny,nz,nc)
    """
    mdata = fmask(data, mask).transpose()

    X = np.atleast_2d(X)
    if X.shape[0] == 1:
        X = X.T

    if isinstance(add_const, bool) and add_const == True:
        basemat = np.atleast_2d(np.ones(np.min(mdata.shape))).T
        nb = 1
    elif isinstance(add_const, np.ndarray) and add_const.shape[1] == X.shape[0]:
        basemat = np.atleast_2d(add_const)
        nb = basemat.shape[0]
    else:
        nb = 0

    if nb != 0:
        X = np.hstack([X, basemat.T])

    tmpbetas = np.linalg.lstsq(X, mdata, rcond=None)[0].transpose()

    if nb != 0:
        tmpbetas = tmpbetas[:, :-nb]

    if betas5d is not None:
        assert Ne5d is not None
        betas5d[:, :, :, Ne5d, :] = unmask(tmpbetas, mask)

    else:
        out = unmask(tmpbetas, mask)
        return out


get_coeffs = _get_coeffs_uni
try:
    if int(os.environ["OMP_NUM_THREADS"]) > 1:
        get_coeffs = _get_coeffs_para
except BaseException:
    pass
