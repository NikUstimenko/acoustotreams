"""ACOUTREAMS: A Python package for acoustic scattering based on the T-matrix method.

.. currentmodule:: acoutreams

Classes
=======

The top-level classes and functions allow a high-level access to the functionality.

Basis sets
----------

   ScalarCylindricalWaveBasis
   ScalarPlaneWaveBasisByUnitVector
   ScalarPlaneWaveBasisByComp
   ScalarSphericalWaveBasis

Matrices and Arrays
-------------------

   AcousticPhysicsArray
   AcousticSMatrix
   AcousticSMatrices
   AcousticTMatrix
   AcousticTMatrixC

Other
-----

   Lattice
   AcousticMaterial

Functions
=========

   vfield
   pfield
   ffamplitude
   expand
   expandlattice
   permute
   plane_wave
   rotate
   translate

"""

_version__ = "0.1.25"

from scipy.special import (  # noqa: F401
    hankel1,
    hankel2,
    jv,
    spherical_jn,
    spherical_yn,
    yv,
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

from treams.misc import(  # noqa: F401
    pickmodes,
    wave_vec_z,
    firstbrillouin1d,
    firstbrillouin2d,
    firstbrillouin3d,
)

from treams._lattice import *   # noqa: F401

from acoutreams._wavesacoustics import *  # noqa: F401

from acoutreams._materialacoustics import AcousticMaterial  # noqa: F401

from acoutreams._smatrixacoustics import (  # noqa: F401
    AcousticSMatrices,
    AcousticSMatrix,
    poynting_avg_z,
)

from acoutreams.scw import *

from acoutreams._coreacoustics import (  # noqa: F401
    ScalarCylindricalWaveBasis,
    AcousticPhysicsArray,
    ScalarPlaneWaveBasisByComp,
    ScalarPlaneWaveBasisByUnitVector,
    ScalarSphericalWaveBasis,
)

from acoutreams._tmatrixacoustics import (  # noqa: F401
    AcousticTMatrix,
    AcousticTMatrixC,
    cylindrical_wave_scalar,
    plane_wave_scalar,
    plane_wave_angle_scalar,
    spherical_wave_scalar,
)

from acoutreams.coeffsacoustics import *  # noqa: F401 

from acoutreams.spw import *  # noqa: F401

from acoutreams._operatorsacoustics import (  # noqa: F401
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

from acoutreams.ssw import *
