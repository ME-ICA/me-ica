#!/usr/bin/env python

import numpy as np
import meica.utils.utils as utils
import pytest


Ne = 3
uncat_shape = (10, 11, 12, Ne, 20)
cat_shape   = (10, 11, 12 * Ne, 20)

np.random.seed(1234)
data_5d = np.random.rand(*uncat_shape)
data_4d = np.random.rand(*cat_shape)
data_2d = np.random.rand(10, 10)


def test_uncat_cat2echos():
    cat = utils.uncat2echos(data_5d, Ne)
    assert cat.shape == cat_shape

    uncat = utils.cat2echos(cat, Ne)
    assert uncat.shape == uncat_shape
    assert np.allclose(uncat, data_5d)

    no_nt = data_5d[:, :, :, :, 0]
    cat = utils.uncat2echos(no_nt, Ne)
    uncat = utils.cat2echos(cat, Ne)

    assert uncat.shape == no_nt.shape + (1,)
    assert np.allclose(np.squeeze(uncat), no_nt)


def test_makemask():
    mask = utils.makemask(data_5d)
    assert mask.shape == uncat_shape[:3]


def test_makeadmask():
    minmask = utils.makeadmask(data_5d, minimum=True)
    assert minmask.shape == uncat_shape[:3]

    minmask = utils.makeadmask(data_5d, minimum=False)
    assert minmask.shape == uncat_shape[:3]

    minmask, masksum = utils.makeadmask(data_5d, minimum=False, getsum=True)
    assert minmask.shape == uncat_shape[:3]


def test_fmask_unmask():
    f_shape = (np.product(cat_shape[:3]),) + cat_shape[3:]
    mask = np.ones(cat_shape[:3])

    fmasked = utils.fmask(data_4d, mask)
    assert fmasked.shape == f_shape

    umasked = utils.unmask(fmasked, mask)
    assert umasked.shape == cat_shape
    assert np.allclose(data_4d, umasked)

    umasked = utils.unmask(fmasked[:, 0], mask)
    assert umasked.shape == cat_shape[:-1]


def test_moments():
    out = utils.moments(data_2d)

    for f in out:
        assert isinstance(f, float)


def test_gaussian():
    h, cx, cy, wx, wy = 1, 1, 1, 1, 1
    out_func = utils.gaussian(h, cx, cy, wx, wy)

    assert out_func(1, 1) == 1


def test_fitgaussian():
    params = utils.fitgaussian(data_2d)

    assert params.shape == (5,)


def test_dice():
    a1, a2 = np.arange(Ne), np.arange(Ne)

    d = utils.dice(a1, a2)
    assert isinstance(d, float)

    d = utils.dice(np.zeros(Ne), np.zeros(Ne))
    assert d == 0

    with pytest.raises(ValueError):
        utils.dice(a1, np.arange(Ne + 1))


def test_andb():
    a1, a2 = np.arange(Ne), np.arange(Ne)

    out = utils.andb([a1, a2])
    assert out.shape == (Ne,)

    with pytest.raises(ValueError):
        utils.andb([a1, np.arange(Ne + 1)])
