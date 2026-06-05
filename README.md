[![PyPI](https://img.shields.io/pypi/v/acoustotreams)](https://pypi.org/project/acoustotreams)
![License](https://img.shields.io/github/license/NikUstimenko/acoustotreams)
[![Documentation Status](https://img.shields.io/badge/docs-online-brightgreen)](https://NikUstimenko.github.io/acoustotreams)
![doctests](https://github.com/NikUstimenko/acoustotreams/actions/workflows/doctests.yml/badge.svg?branch=main)
![tests](https://github.com/NikUstimenko/acoustotreams/actions/workflows/tests.yml/badge.svg?branch=main)
[![coverage](https://img.shields.io/endpoint?url=https://NikUstimenko.github.io/acoustotreams/endpoint.json)](https://htmlpreview.github.io/?https://github.com/NikUstimenko/acoustotreams/blob/htmlcov/index.html)

# acoustotreams

The package `acoustotreams` adopts the framework of the `treams` package for pressure-wave scattering in arrangements of multiple acoustic scatterers, with and without peridic boundary conditions, using the T-matrix method.

## Installation

### Installation using pip

To install the package with pip, use

```sh
pip install acoustotreams
```
Preliminarily, you have to also install original `treams` as well as `numpy<2` and `scipy<1.11`
```sh
pip install treams
```

## Documentation

The documentation can be found at https://NikUstimenko.github.io/acoustotreams.

## Publications

When using this code please cite:

The following publications document the developments and methods for different parts of the code:

* [N. Ustimenko, I. Fernandez-Corbaton, and C. Rockstuhl, Singular value decomposition to describe bound states in the continuum in periodic metasurfaces, arXiv 2602.15741 (2026).](https://arxiv.org/abs/2602.15741)
* [O. Demeulenaere, N. Ustimenko, A. G. Athanassiadis,  L. Gulati, C. Rockstuhl, and P. Fischer, Ultrasonic metamaterial at MHz frequencies using microstructured glass, arXiv 2512.20506 (2026).](https://arxiv.org/abs/2512.20506)
* [N. Ustimenko, A. B. Evlyukhin, V. Kyrimi, A. V. Kildishev, and C. Rockstuhl, Lattice-induced sound trapping in biperiodic metasurfaces of acoustic resonators, Phys. Rev. Res. 8, 013074 (2026).](https://doi.org/10.1103/wnmk-zhrb)
* [N. Ustimenko, C. Rockstuhl, and A. V. Kildishev, Optimal multipole center for subwavelength acoustic scatterers, Appl. Phys. Lett. 126, 142201 (2025).](https://doi.org/10.1063/5.0257760)

## Features

* [x] T-matrix calculations using a spherical or cylindrical wave basis set
* [x] Scattering from clusters of particles
* [x] Scattering from particles and clusters arranged in 3d-, 2d-, and 1d-lattices
* [x] Calculation of sound propagation in stratified media
* [x] Band calculation in crystal structures
