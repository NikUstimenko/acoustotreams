r"""Special (mathematical) functions.

.. currentmodule:: acoustotreams.special

Special mathematical functions used in :mod:`acoustotreams`. Some functions are imported from
:py:mod:`scipy.special` and :py:mod:`treams.special`.

Available functions
===================

Spherical waves
---------------

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
-----------------

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
-----------

.. autosummary::
   :toctree: generated/

   spw_Psi
   vpw_L


Those functions are just imported from Scipy. So, one only needs to import this
package within acoustotreams.

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

The following functions are just imported from treams. So, one only needs to import this
package within treams.

Bessel and Hankel functions, with their spherical counterparts, derivatives
---------------------------------------------------------------------------

:py:data:`~treams.special.hankel1_d`
:py:data:`~treams.special.hankel2_d`
:py:data:`~treams.special.jv_d`
:py:data:`~treams.special.yv_d`
:py:data:`~treams.special.spherical_hankel1`
:py:data:`~treams.special.spherical_hankel2`
:py:data:`~treams.special.spherical_hankel1_d`
:py:data:`~treams.special.spherical_hankel2_d`

Those functions just wrap Scipy functions with special optional arguments to be able to
analogously access them like their non-spherical counterparts:

:py:data:`~treams.special.spherical_jn_d`
:py:data:`~treams.special.spherical_yn_d`


Scipy functions with enhanced domain
------------------------------------

:py:data:`~treams.special.sph_harm`
:py:data:`~treams.special.lpmv`


Integrals for the Ewald summation
---------------------------------

:py:data:`~treams.special.incgamma`
:py:data:`~treams.special.intkambe`


Wigner d- and Wigner D-matrix elements
--------------------------------------

:py:data:`~treams.special.wignersmalld`
:py:data:`~treams.special.wignerd`


Wigner 3j-symbols
-----------------

:py:data:`~treams.special.wigner3j`

Vector spherical harmonics
-----------------
   
:py:data:`~treams.special.pi_fun`
:py:data:`~treams.special.tau_fun`
:py:data:`~treams.special.vsh_X`
:py:data:`~treams.special.vsh_Y`
:py:data:`~treams.special.vsh_Z`


Coordinate system transformations
---------------------------------

:py:data:`~treams.special.car2cyl`
:py:data:`~treams.special.car2sph`
:py:data:`~treams.special.cyl2car`
:py:data:`~treams.special.cyl2sph`
:py:data:`~treams.special.sph2car`
:py:data:`~treams.special.sph2cyl`
:py:data:`~treams.special.vcar2cyl`
:py:data:`~treams.special.vcar2sph`
:py:data:`~treams.special.vcyl2car`
:py:data:`~treams.special.vcyl2sph`
:py:data:`~treams.special.vsph2car`
:py:data:`~treams.special.vsph2cyl`
:py:data:`~treams.special.car2pol`
:py:data:`~treams.special.pol2car`
:py:data:`~treams.special.vcar2pol`
:py:data:`~treams.special.vpol2car`

"""

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

from acoustotreams.special._wavesacoustics import *  # noqa: F401