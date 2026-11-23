[![PyPI](https://img.shields.io/pypi/v/acoustotreams)](https://pypi.org/project/acoustotreams)
![License](https://img.shields.io/github/license/NikUstimenko/acoustotreams)
[![Documentation Status](https://img.shields.io/badge/docs-online-brightgreen)](https://NikUstimenko.github.io/acoustotreams)
![doctests](https://github.com/NikUstimenko/acoustotreams/actions/workflows/doctests.yml/badge.svg?branch=main)
![tests](https://github.com/NikUstimenko/acoustotreams/actions/workflows/tests.yml/badge.svg?branch=main)
[![coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fraw.githubusercontent.com%2FNikUstimenko%2Facoustotreams%2Fhtmlcov%2Fendpoint.json)](https://github.com/NikUstimenko/acoustotreams/tree/htmlcov)

# acoustotreams

The package `acoustotreams` adopts the framework of the `treams` package for pressure-wave scattering in arrangements of multiple acoustic scatterers, with and without periodic boundary conditions, using the T-matrix method. 

## Installation

### Installation using pip

`acoustotreams` supports Python 3.9–3.13. We recommend installing it in a [virtual environment](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/) or a [conda environment]( https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Create and activate one with a supported Python version.

Inside an activated environment, install the package with

```sh
pip install acoustotreams
```

or

```sh
python -m pip install acoustotreams
```

Required dependencies (`treams`, `numpy`, `scipy`) are installed automatically. To update later, run `python -m pip install --upgrade acoustotreams`.

> **Note:** Python 3.14 and newer are not supported yet. With these versions, pip falls back to an outdated release and the installation fails with a compiler error. If your default Python is 3.14 or newer, create the environment with an older version.

> **Note (Windows):** On Windows, use Python 3.9–3.12. `treams` does not work with NumPy 2 on Windows yet, and NumPy 1.x is not available for Python 3.13.

If you are not using a virtual environment, invoke pip through a supported interpreter directly, for example Python 3.12.

On macOS or Linux:

```sh
python3.12 -m pip install acoustotreams
```

On Windows:

```bat
py -3.12 -m pip install acoustotreams
```

To verify the installation, run

```sh
python -c "import acoustotreams; from importlib.metadata import version; print('acoustotreams', version('acoustotreams'))"
```

Outside a virtual environment, replace `python` with `python3.13` (macOS / Linux) or `py -3.13` (Windows). The command should print the latest version listed on [PyPI](https://pypi.org/project/acoustotreams/).

#### Other platforms

Pre-built `treams` packages are available for Linux (x86_64), macOS (Apple Silicon) and Windows (64-bit). There are none for other platforms, such as Intel Macs or ARM-based Linux, and pip cannot install a working version of `treams` there automatically. Build `treams` from source first (requires git and a C compiler), then install `acoustotreams`.


## Documentation

The documentation can be found at https://NikUstimenko.github.io/acoustotreams.

## Publications

When using this code, please cite:

[N. Ustimenko, and C. Rockstuhl, acoustotreams -- A Python package for acoustic-wave scattering based on the T-matrix method, arXiv 2606.22573 (2026).](https://arxiv.org/abs/2606.22573)

The following publications document the developments and methods for different parts of the code:

* [N. Ustimenko, I. Fernandez-Corbaton, and C. Rockstuhl, Singular value decomposition to describe bound states in the continuum in periodic metasurfaces, Phys. Rev. B 113, 235418 (2026).](https://doi.org/10.1103/gx73-kyj9)
* [O. Demeulenaere, N. Ustimenko, A. G. Athanassiadis,  L. Gulati, C. Rockstuhl, and P. Fischer, Ultrasonic metamaterial at MHz frequencies using microstructured glass, Science Advances 12, eaef0666 (2026).](https://doi.org/10.1126/sciadv.aef0666)
* [N. Ustimenko, A. B. Evlyukhin, V. Kyrimi, A. V. Kildishev, and C. Rockstuhl, Lattice-induced sound trapping in biperiodic metasurfaces of acoustic resonators, Phys. Rev. Res. 8, 013074 (2026).](https://doi.org/10.1103/wnmk-zhrb)
* [N. Ustimenko, C. Rockstuhl, and A. V. Kildishev, Optimal multipole center for subwavelength acoustic scatterers, Appl. Phys. Lett. 126, 142201 (2025).](https://doi.org/10.1063/5.0257760)

## Features

* [x] T-matrix calculations using a spherical or cylindrical wave basis set
* [x] Scattering from clusters of particles
* [x] Scattering from particles and clusters arranged in 3d-, 2d-, and 1d-lattices
* [x] Calculation of sound propagation in stratified media
* [x] Band calculation in crystal structures
