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

Spherical waves and translation coefficients

.. autosummary::
   :toctree:

   ssw_Psi
   ssw_rPsi

Cylindrical waves

.. autosummary::
   :toctree:

Plane waves

.. autosummary::
   :toctree:

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

.. autosummary::
   :toctree: generated/
   :nosignatures:

   treams.special.spherical_jn_d
   treams.special.spherical_yn_d
   treams.special.sph_harm
   treams.special.lpmv
   treams.special.incgamma
   treams.special.intkambe
   treams.special.wignersmalld
   treams.special.wignerd
   treams.special.wigner3j
   treams.special.pi_fun
   treams.special.tau_fun
   treams.special.vsh_X
   treams.special.vsh_Y
   treams.special.vsh_Z
   treams.special.car2cyl
   treams.special.car2sph
   treams.special.cyl2car
   treams.special.cyl2sph
   treams.special.sph2car
   treams.special.sph2cyl
   treams.special.vcar2cyl
   treams.special.vcar2sph
   treams.special.vcyl2car
   treams.special.vcyl2sph
   treams.special.vsph2car
   treams.special.vsph2cyl
   treams.special.car2pol
   treams.special.pol2car
   treams.special.vcar2pol
   treams.special.vpol2car
   treams.misc.wave_vec_z
   treams.misc.firstbrillouin1d
   treams.misc.firstbrillouin2d
   treams.misc.firstbrillouin3d
   treams.lattice



"""

_version__ = "0.1.17"

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

from acoustotreams.scw import *

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

from acoustotreams.coeffsacoustics import *  # noqa: F401 

from acoustotreams.spw import *  # noqa: F401

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

from acoustotreams.ssw import *
