"""ACOUSTOTREAMS: A Python package for acoustic wave scattering based on the T-matrix method.

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

Other
-----

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
   vamplitudeff
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

Spherical waves

.. autosummary::
   :toctree: generated/

   ssw_Psi
   ssw_psi
   vsw_L
   vsw_l
   ssw_rPsi
   vsw_rL
   tl_ssw
   tl_ssw_r

Cylindrical waves

.. autosummary::
   :toctree: generated/

   scw_Psi
   scw_psi
   vcw_L
   vcw_l
   scw_rPsi
   vcw_rL
   tl_scw
   tl_scw_r


Plane waves

.. autosummary::
   :toctree: generated/

   spw_Psi
   vpw_L

Functions imported from SciPy
-----------------------------

+------------------------------------------------------------+-------------------------+
| :py:data:`~scipy.special.hankel1`\(v, z[, out])            | Hankel function of the  |
|                                                            | first kind.             |
+------------------------------------------------------------+-------------------------+
| :py:data:`~scipy.special.hankel2`\(v, z[, out])            | Hankel function of the  |
|                                                            | second kind.            |
+------------------------------------------------------------+-------------------------+
| :py:data:`~scipy.special.jv`\(v, z[, out])                 | Bessel function of the  |
|                                                            | first kind of real      |
|                                                            | order and complex       |
|                                                            | argument.               |
+------------------------------------------------------------+-------------------------+
| :py:data:`~scipy.special.yv`\(v, z[, out])                 | Bessel function of the  |
|                                                            | second kind of real     |
|                                                            | order and complex       |
|                                                            | argument.               |
+------------------------------------------------------------+-------------------------+
| | :py:func:`spherical_jn <scipy.special.spherical_jn>`\(n, | Spherical Bessel        |
|   z[, derivative])                                         | function of the first   |
|                                                            | kind or its derivative. |
+------------------------------------------------------------+-------------------------+
| | :py:func:`spherical_yn <scipy.special.spherical_yn>`\(n, | Spherical Bessel        |
|   z[, derivative])                                         | function of the second  |
|                                                            | kind or its derivative. |
+------------------------------------------------------------+-------------------------+

Functions imported from treams.special, treams.misc, and treams.lattice
-----------------------------------------------------------------------

+------------------------------------------------------------+-----------------------------+
| :py:data:`~treams.special.spherical_jn_d`\(n, z)           | Derivative of the spherical | 
|                                                            | Bessel function of the      |
|                                                            | first kind.                 |   
+------------------------------------------------------------+-----------------------------+
| :py:data:`~treams.special.spherical_yn_d`\(n, z)           | Derivative of the spherical | 
|                                                            | Bessel function of the      |
|                                                            | second kind.                |   
+------------------------------------------------------------+-----------------------------+

"""

_version__ = "0.1.38"

from scipy.special import (  # noqa: F401
    hankel1,
    hankel2,
    jv,
    yv,
    spherical_jn,
    spherical_yn,
)

from treams.special import(   # noqa: F401
    spherical_jn_d,
    spherical_yn_d,
    sph_harm,
    lpmv,
    incgamma,
    intkambe,
    wignersmalld,
    wignerd,
    wigner3j,
    pi_fun,
    tau_fun,
    vsh_X,
    vsh_Y,
    vsh_Z,
    car2cyl,
    car2sph,
    cyl2car,
    cyl2sph,
    sph2car,
    sph2cyl,
    vcar2cyl,
    vcar2sph,
    vcyl2car,
    vcyl2sph,
    vsph2car,
    vsph2cyl,
    car2pol,
    pol2car,
    vcar2pol,
    vpol2car,
    )

from treams.misc import(  # noqa: F401,
    wave_vec_z,
    firstbrillouin1d,
    firstbrillouin2d,
    firstbrillouin3d,
)

from treams._lattice import *   # noqa: F401

from acoustotreams._wavesacoustics import *  # noqa: F401

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
    VAmplitudeFF,
    Expand,
    ExpandLattice,
    Permute,
    Rotate,
    Translate,
    vfield,
    pfield,
    pamplitudeff,
    vamplitudeff,
    expand,
    expandlattice,
    permute,
    rotate,
    translate,
)

import acoustotreams.ssw # noqa: F401
