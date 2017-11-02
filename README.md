# ME-ICA

Core code for Multi-Echo Independent Components Analysis (ME-ICA).

[![Build Status](https://travis-ci.org/emdupre/me-ica.svg?branch=shotgun)](https://travis-ci.org/emdupre/me-ica) [![CircleCI](https://circleci.com/gh/emdupre/me-ica/tree/py3.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/emdupre/me-ica)
[![Documentation Status](https://readthedocs.org/projects/me-ica/badge/?version=latest)](http://me-ica.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/emdupre/me-ica/branch/shotgun/graph/badge.svg)](https://codecov.io/gh/emdupre/me-ica) [![License](https://img.shields.io/badge/License-LGPL%202.0-blue.svg)](https://opensource.org/licenses/LGPL-2.1)

## About

`meica` is a pipeline for preprocessing multi-echo datasets and applying multi-echo ICA to the spatially concatenated echoes. It includes the following steps:

1. Calculates motion parameters based on images with highest contrast (usually the first echo)
2. Applies motion correction and T2*-weighted co-registration parameters
3. Applies standard EPI preprocessing (slice-time correction, etc.)
4. Computes PCA and ICA in conjunction with TE-dependence analysis

For more information, please see the documentation:

https://me-ica.readthedocs.io/

## Contributors

| [<img src="https://avatars.githubusercontent.com/emdupre?s=100" width="100" alt="Elizabeth DuPre" /><br /><sub>Elizabeth DuPre</sub>](http://emdupre.me)<br />[💻 📖 ⚠️](https://github.com/emdupre/me-ica/commits?author=emdupre) | [<img src="https://avatars.githubusercontent.com/prantikk?s=100" width="100" alt="Prantik Kundu" /><br /><sub>Prantik Kundu</sub>](http://www.mountsinai.org/profiles/prantik-kundu)<br />[💻](https://github.com/emdupre/me-ica/commits?author=prantikk) |[<img src="https://avatars.githubusercontent.com/rmarkello?s=100" width="100" alt="Ross Markello" /><br /><sub>Ross Markello</sub>](http://rossmarkello.me/)<br />[📖 ⚠️](https://github.com/emdupre/me-ica/commits?author=rmarkello) |[<img src="https://avatars.githubusercontent.com/tsalo?s=100" width="100" alt="Taylor Salo" /><br /><sub>Taylor Salo</sub>](https://tsalo.github.io/)<br />[💻 📖](https://github.com/emdupre/me-ica/commits?author=tsalo) | [<img src="https://avatars.githubusercontent.com/KirstieJane?s=100" width="100" alt="Kirstie Whitaker" /><br /><sub>Kirstie Whitaker</sub>](http://whitakerlab.github.io)<br />[📖](https://github.com/emdupre/me-ica/commits?author=KirstieJane)
| :---: | :---: | :---: | :---: | :---: |


This project follows the [all-contributors](https://github.com/kentcdodds/all-contributors) specification.

## Citations
When using ME-ICA, please include the following citations:

> Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J., Saad, Z.S., Bandettini, P.A. & Bullmore, E.T. (2013). Integrated strategy for improving functional connectivity mapping using multiecho fMRI. Proceedings of the National Academy of Sciences, 110, 16187–16192.



> Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A. (2011). Differentiating BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage, 60, 1759-1770.
