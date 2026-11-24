r"""ACOUSTOTREAMS: A Python package for acoustic wave scattering based on the T-matrix method.

.. currentmodule:: acoustotreams

Classes
=======

The top-level classes and functions allow a high-level access to the functionality.

Basis sets
----------

.. autosummary::
   :toctree: generated/

   ScalarCylindricalWaveBasis
   ScalarPlaneWaveBasisByUnitVector
   ScalarPlaneWaveBasisByComp
   ScalarSphericalWaveBasis

Matrices and Arrays
-------------------

.. autosummary::
   :toctree: generated/

   AcousticsArray
   AcousticSMatrix
   AcousticSMatrices
   AcousticTMatrix
   AcousticTMatrixC

Other classes
-------------

.. autosummary::
   :toctree: generated/

   AcousticMaterial

Functions
=========

Operators
---------

.. autosummary::
   :toctree: generated/

   pfield
   vfield
   pamplitudeff
   expand
   expandlattice
   permute
   rotate
   translate

Scalar wave functions
---------------------

.. autosummary::
   :toctree: generated/

   cylindrical_wave_scalar
   plane_wave_scalar
   plane_wave_angle_scalar
   spherical_wave_scalar

Other functions
---------------

.. autosummary::
   :toctree: generated/

   poynting_avg_z

Functions and classes imported from treams
------------------------------------------

Special functions, including those imported from SciPy and treams, are available in
:mod:`acoustotreams.special`.

+--------------------------------------------------------------+-------------------------------+
| :py:func:`~treams.misc.wave_vec_z`\(kx, ky, k)               | Z component of the wave       |
|                                                              | vector with positive          |
|                                                              | imaginary part.               |
+--------------------------------------------------------------+-------------------------------+
| :py:func:`~treams.misc.firstbrillouin1d`\(kpar, b)           | Map wave vector to the first  |
|                                                              | Brillouin zone in 1D.         |
+--------------------------------------------------------------+-------------------------------+
| :py:func:`~treams.misc.firstbrillouin2d`\(kpar, b, n=2)      | Map wave vector to the first  |
|                                                              | Brillouin zone in 2D.         |
+--------------------------------------------------------------+-------------------------------+
| :py:func:`~treams.misc.firstbrillouin3d`\(kpar, b, n=2)      | Map wave vector to the first  |
|                                                              | Brillouin zone in 3D.         |
+--------------------------------------------------------------+-------------------------------+
| :py:class:`~treams.Lattice`\(arr, alignment=None)            | Lattice definition in one,    |
|                                                              | two, or three dimensions.     |
+--------------------------------------------------------------+-------------------------------+
| :py:class:`WaveVector`\(seq=(), alignment=None)      | Wave vector with up to three  |
|                                                              | components and an alignment.  |
+--------------------------------------------------------------+-------------------------------+

"""

__version__ = "0.2.53"

import sys

import numpy as np

if sys.platform == "win32" and int(np.__version__.split(".")[0]) >= 2:
    raise ImportError(
        "treams does not support NumPy 2 on Windows yet. Use Python 3.12 or "
        "older, where acoustotreams installs NumPy 1.x automatically."
    )

from treams.misc import(  # noqa: F401,
    wave_vec_z,
    firstbrillouin1d,
    firstbrillouin2d,
    firstbrillouin3d,
)

from treams._lattice import Lattice, WaveVector  # noqa: F401

from acoustotreams._materialacoustics import AcousticMaterial  # noqa: F401

from acoustotreams._smatrixacoustics import (  # noqa: F401
    AcousticSMatrices,
    AcousticSMatrix,
    poynting_avg_z,
)

import acoustotreams.scw  # noqa: F401

from acoustotreams._coreacoustics import (  # noqa: F401
    ScalarCylindricalWaveBasis,
    AcousticsArray,
    ScalarPlaneWaveBasisByComp,
    ScalarPlaneWaveBasisByUnitVector,
    ScalarSphericalWaveBasis,
)

from acoustotreams._tmatrixacoustics import (  # noqa: F401
    AcousticTMatrix,
    AcousticTMatrixC,
    cylindrical_wave_scalar,
    plane_wave_scalar,
    plane_wave_angle_scalar,
    spherical_wave_scalar,
)

import acoustotreams.coeffs  # noqa: F401 

import acoustotreams.spw  # noqa: F401

from acoustotreams._operatorsacoustics import (  # noqa: F401
    PField,
    VField,
    PAmplitudeFF,
    Expand,
    ExpandLattice,
    Permute,
    Rotate,
    Translate,
    vfield,
    pfield,
    pamplitudeff,
    expand,
    expandlattice,
    permute,
    rotate,
    translate,
)

import acoustotreams.ssw # noqa: F401
