import warnings

import numpy as np
import sklearn.decomposition._fastica
from scipy.stats import spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.extmath import check_random_state

_sym_decorrelation = sklearn.decomposition._fastica._sym_decorrelation

INITIAL_TOL_OK = 0.1
CHECK_STEPS_OK = 150
INITIAL_CHECK_FAC = 5
TREND_SIG = 0.2


def init_w(nc):
    random_state = check_random_state(None)
    w_init = np.asarray(random_state.normal(size=(nc, nc)), dtype=float)
    return w_init


def steps_ok(steps_lim, acc_lim=INITIAL_TOL_OK):
    all_check = np.array(steps_lim) < acc_lim
    if np.sum(all_check.astype(int)) == len(steps_lim):
        return True
    else:
        return False


def is_converging(lim_list, initially=False):
    lim_list = np.array(lim_list)
    if initially:
        if steps_ok(lim_list[-(CHECK_STEPS_OK // INITIAL_CHECK_FAC) :]):
            return True
        else:
            return False
    else:
        test_sam = lim_list[-CHECK_STEPS_OK:]
        _res = spearmanr(np.arange(len(test_sam)), test_sam)
        if _res.statistic < 0 and _res.pvalue < TREND_SIG:
            print(f"Convergence significance: {_res.pvalue}")
            return True
        else:
            print(f"Detected non-ideal convergence: {_res}")
            return False


def _ica_par_dyn(X, tol, g, fun_args, max_iter, w_init, n_attempts=5):
    """Dynamic Parallel FastICA.

    Aborts and retarts given parameters of convergence

    Used internally by FastICA --main loop.
    To be monkey patched with new features
      - reporting values
      - restarting
    """
    done = False
    min_restart_lim = 50 * tol
    W = _sym_decorrelation(w_init)
    del w_init
    p_ = float(X.shape[1])
    ii = -1
    for att in range(n_attempts):
        if done:
            break
        initial_conv = False
        lim_his = []
        if att != 0:
            w_init = init_w(nc=X.shape[0])
            W = _sym_decorrelation(w_init)
            del w_init
            print(f"Attempt {att}")
        for ii in range(max_iter):
            restart = False
            gwtx, g_wtx = g(np.dot(W, X), fun_args)
            W1 = _sym_decorrelation(np.dot(gwtx, X.T) / p_ - g_wtx[:, np.newaxis] * W)
            del gwtx, g_wtx
            # builtin max, abs are faster than numpy counter parts.
            # np.einsum allows having the lowest memory footprint.
            # It is faster than np.diag(np.dot(W1, W.T)).
            lim = max(abs(abs(np.einsum("ij,ij->i", W1, W)) - 1))
            lim_his.append(lim)
            print("Step %i: %.9f" % (ii, lim))
            W = W1

            if lim < tol:
                done = True
                break

            if not done and len(lim_his) > CHECK_STEPS_OK and ii % 5 == 0:
                if ii > max_iter:
                    restart = True
                elif not initial_conv:
                    if is_converging(lim_his, initially=True):
                        initial_conv = True
                    else:
                        restart = True
                else:
                    if not is_converging(lim_his) and lim > min_restart_lim:
                        restart = True

                if restart:
                    warnings.warn("Restarting optimization.")
                    break
                else:
                    continue

    if not done:
        warnings.warn(
            "FastICA did not converge. Consider increasing "
            "tolerance or the maximum number of iterations.",
            ConvergenceWarning,
        )

    np.savetxt(
        "convergence_history.1D",
        np.array(lim_his),
        fmt="%0.7f",
        delimiter=" ",
        header=f"# Converged: {done}   Steps: {ii + 1}",
    )

    return (
        W,
        lim,
    )  # returning lim is hack to return final tol,
    #     instead of num_iter which FastICA.transform stores


sklearn.decomposition._fastica._ica_par = _ica_par_dyn
FastICA = sklearn.decomposition._fastica.FastICA
