[![PyPI](https://img.shields.io/pypi/v/acoustotreams)](https://pypi.org/project/acoustotreams)
![License](https://img.shields.io/github/license/NikUstimenko/acoustotreams)
[![Documentation Status](https://img.shields.io/badge/docs-online-brightgreen)](https://NikUstimenko.github.io/acoustotreams)
![doctests](https://github.com/NikUstimenko/acoustotreams/actions/workflows/doctests.yml/badge.svg?branch=main)
![tests](https://github.com/NikUstimenko/acoustotreams/actions/workflows/tests.yml/badge.svg?branch=main)
[![coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fraw.githubusercontent.com%2FNikUstimenko%2Facoustotreams%2Fhtmlcov%2Fendpoint.json)](https://github.com/NikUstimenko/acoustotreams/tree/htmlcov)

# acoustotreams

The package `acoustotreams` adopts the framework of the `treams` package for pressure-wave scattering in arrangements of multiple acoustic scatterers, with and without peridic boundary conditions, using the T-matrix method. 

## Installation

### Installation using pip

`acoustotreams` requires Python 3.9–3.13. To install the package with pip, use

```sh
pip install acoustotreams
```

Alternatively, on macOS or Linux, invoke pip through Python 3:

```sh
python3 -m pip install acoustotreams
```

This runs pip with the selected Python 3 interpreter. Required dependencies, including `treams`, are installed automatically.

On Windows, use:

```sh
py -m pip install acoustotreams
```

To verify the installation on macOS or Linux, run:

```sh
python3 -c "import acoustotreams; print('acoustotreams imported successfully')"
```

On Windows, replace `python3` with `py` in the verification command.

We recommend installing in a [virtual environment](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/).

## Documentation

The documentation can be found at https://NikUstimenko.github.io/acoustotreams.

## Publications

When using this code, please cite:

[N. Ustimenko, and C. Rockstuhl, acoustotreams -- A Python package for acoustic-wave scattering based on the T-matrix method, arXiv 2606.22573 (2026).](https://arxiv.org/abs/2606.22573)

The following publications document the developments and methods for different parts of the code:

* [N. Ustimenko, I. Fernandez-Corbaton, and C. Rockstuhl, Singular value decomposition to describe bound states in the continuum in periodic metasurfaces, Phys. Rev. B 113, 235418 (2026).](https://doi.org/10.1103/gx73-kyj9)
* [O. Demeulenaere, N. Ustimenko, A. G. Athanassiadis,  L. Gulati, C. Rockstuhl, and P. Fischer, Ultrasonic metamaterial at MHz frequencies using microstructured glass, Science Advances 12, aef0666 (2026).](https://doi.org/10.1126/sciadv.aef0666)
* [N. Ustimenko, A. B. Evlyukhin, V. Kyrimi, A. V. Kildishev, and C. Rockstuhl, Lattice-induced sound trapping in biperiodic metasurfaces of acoustic resonators, Phys. Rev. Res. 8, 013074 (2026).](https://doi.org/10.1103/wnmk-zhrb)
* [N. Ustimenko, C. Rockstuhl, and A. V. Kildishev, Optimal multipole center for subwavelength acoustic scatterers, Appl. Phys. Lett. 126, 142201 (2025).](https://doi.org/10.1063/5.0257760)

## Features

* [x] T-matrix calculations using a spherical or cylindrical wave basis set
* [x] Scattering from clusters of particles
* [x] Scattering from particles and clusters arranged in 3d-, 2d-, and 1d-lattices
* [x] Calculation of sound propagation in stratified media
* [x] Band calculation in crystal structures
