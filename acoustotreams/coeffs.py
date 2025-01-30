"""Scattering coefficients for high-symmetry cases.

Calculate the scattering coefficients for cases where they can be obtained analytically
easily. This is a sphere using spherical waves (Mie coefficients), a
cylinder using cylindrical waves, and an infinitely extended planar
interface (Fresnel coefficients).

.. autosummary::
   :toctree:

   mie_acoustics
   mie_acoustics_cyl
   fresnel_acoustics

"""

import treams.special
import numpy as np
from acoustotreams._materialacoustics import AcousticMaterial

def mie_acoustics(l, x, *materials):
    r"""Mie scattering coefficient of degree l.

    The sphere is defined by its size parameter :math:`k_0 r`, where :math:`r` is the
    radius and :math:`k_0` the wave number in air. 

    Likewise, the material parameters are given from inside to outside. These arrays
    are expected to be exactly one unit larger then the array `x`.

    The result is a complex number relating incident with the scattered modes, which are 
    index in the same way.

    Args:
        l (integer): Degree :math:`l \geq 0`
        x (float, array_like): Size parameters
        rho (float or complex, array_like): Mass density
        c (float or complex, array_like): Longitudinal speed of sound
        c_t (float or complex, array_like): Transverse speed of sound

    Returns:
        complex
    """
    mat_sphere, mat_env = zip(*materials) 
    x_env = x * AcousticMaterial().c / mat_env[1]
    j = treams.special.spherical_jn(l, x_env, derivative=False)
    j_d = treams.special.spherical_jn(l, x_env, derivative=True)
    h = treams.special.spherical_hankel1(l, x_env)
    h_d = treams.special.spherical_hankel1_d(l, x_env)
    if np.abs(mat_sphere[2]) == 0 and np.abs(mat_sphere[1]) == 0 and np.abs(mat_sphere[0]) == np.inf:
        res = -j_d / h_d
    elif np.abs(mat_sphere[2]) == 0 and np.abs(mat_sphere[1]) == 0 and np.abs(mat_sphere[0]) == 0:
        res = -j / h
    elif np.abs(mat_sphere[1]) > 0:
        x_sphere = x * AcousticMaterial().c / mat_sphere[1]
        j1 = treams.special.spherical_jn(l, x_sphere, derivative=False)
        j1_d = treams.special.spherical_jn(l, x_sphere, derivative=True)
        if np.abs(mat_sphere[2]) == 0:
            delta = x_sphere * mat_env[0] / (x_env * mat_sphere[0]) 
            res = (delta * j1_d * j - j1 * j_d) / (j1 * h_d - delta * j1_d * h) 
        elif np.abs(mat_sphere[2]) > 0:
            x_sphere_t = x * AcousticMaterial().c / mat_sphere[2]
            j1t = treams.special.spherical_jn(l, x_sphere_t, derivative=False)
            j1t_d = treams.special.spherical_jn(l, x_sphere_t, derivative=True)
            if isinstance(l, (tuple, list)):
                l = np.array(l)
            d44 = (l * (l + 1) - 0.5 * x_sphere_t**2) * j1 - 2 * x_sphere * j1_d
            d33 = (l * (l + 1) - 0.5 * x_sphere_t**2 - 1) * j1t - x_sphere_t * j1t_d
            d34 = x_sphere * j1_d - j1
            d43 = l * (l + 1) * (x_sphere_t * j1t_d - j1t)
            d24 = x_sphere * j1_d
            d23 = l * (l + 1) * j1t
            zeta = (d44*d33 - d34*d43) / (d24*d33 - d34*d23)
            res = (-(mat_env[0]/mat_sphere[0] * 0.5 * x_sphere_t**2/x_env * j + zeta * j_d) / 
                (mat_env[0]/mat_sphere[0] * 0.5 * x_sphere_t**2/x_env * h + zeta * h_d))
    return res

def mie_acoustics_cyl(kz, m, k0, radii, *materials):
    r"""Scattering coefficient at an infinite cylinder

    The cylinder is defined by its radii.
    Likewise, the material parameters are given from inside to outside. These arrays
    are expected to be exactly one unit larger then the array `radii`.

    The result is a complex number relating incident with the scattered modes, 
    which are index in the same way.

    Args:
        kz (float): Z component of the wave
        m (integer): Order
        k0 (float or complex): Wave number in vacuum
        radii (float, array_like): Size parameters
        rho (float or complex, array_like): Mass density
        c (float or complex, array_like): Longitudinal speed of sound
        c_t (float or complex, array_like): Transverse speed of sound

    Returns:
        complex
    """

    mat_cyl, mat_env = zip(*materials) 
    k_env = k0 * AcousticMaterial().c / mat_env[1] + 0j
    krho_env = np.sqrt(k_env * k_env - kz * kz)
    x_env = krho_env * radii
    j = treams.special.jv(m, x_env)
    j_d = treams.special.jv_d(m, x_env)
    h = treams.special.hankel1(m, x_env)
    h_d = treams.special.hankel1_d(m, x_env)
    
    if np.abs(mat_cyl[2]) == 0 and np.abs(mat_cyl[1]) == 0 and np.abs(mat_cyl[0]) == np.inf:
        res = -j_d / h_d
    elif np.abs(mat_cyl[2]) == 0 and np.abs(mat_cyl[1]) == 0 and np.abs(mat_cyl[0]) == 0:
        res = -j / h
    if np.abs(mat_cyl[1]) > 0:
        k_cyl = k0 * AcousticMaterial().c / mat_cyl[1] + 0j
        krho_cyl = np.sqrt(k_cyl * k_cyl - kz * kz)
        x_cyl = krho_cyl * radii
        j1 = treams.special.jv(m, x_cyl)
        j1_d = treams.special.jv_d(m, x_cyl)
        if np.abs(mat_cyl[2]) == 0:
            delta = x_cyl * mat_env[0] / (x_env * mat_cyl[0]) 
            res = (delta * j1_d * j - j1 * j_d) / (j1 * h_d - delta * j1_d * h)
        elif np.abs(mat_cyl[2]) > 0:
            k_cyl_t = k0 * AcousticMaterial().c / mat_cyl[2] + 0j
            krho_cyl_t = np.sqrt(k_cyl_t * k_cyl_t - kz * kz) 
            x_cyl_t = krho_cyl_t * radii
            j1t = treams.special.jv(m, x_cyl_t)
            j1t_d = treams.special.jv_d(m, x_cyl_t)
            j1t_dd = (-j1t_d + (x_cyl_t**2 + m**2) * j1t) / x_cyl_t
            j1_dd = (-j1_d + (x_cyl**2 + m**2) * j1) / x_cyl
            matrix = np.zeros((3, 3), complex)
            rhs = np.zeros(3, complex)
            lam_0 = mat_env[0] * mat_env[1]**2
            mu_1 = mat_cyl[0] * mat_cyl[2]**2
            lam_1 = mat_cyl[0] * mat_cyl[1]**2 - 2 * mu_1
            matrix[0][1] = 2j * kz * krho_cyl/k_cyl * j1_d
            matrix[0][2] = (krho_cyl_t**2 - kz**2)/k_cyl_t * j1t_d
            matrix[1][0] = k_env * lam_0 * h
            matrix[1][1] = 2 * mu_1 * krho_cyl**2/k_cyl * j1_dd - k_cyl * lam_1 * j1
            matrix[1][2] = 2j * mu_1 * kz * krho_cyl_t/k_cyl_t * j1t_dd
            rhs[1] = -np.power(-1j, m) * k_env * lam_0 * j
            matrix[2][0] = -krho_env/k_env * h_d
            matrix[2][1] = krho_cyl/k_cyl * j1_d
            matrix[2][2] = 1j * kz/k_cyl_t * j1t_d
            rhs[2] = np.power(-1j, m) * krho_env/k_env * j_d
            res = np.linalg.solve(matrix, rhs)
            res = res[0] * np.power(1j, m)
    return res

def fresnel_acoustics(kzs, rhos):
    r"""Fresnel coefficients for a planar interface.

    The first two dimensions index the two media for the above
    and below the S-matrix, the second two dimensions are added 
    to meet the treams convention.

    The result is an array relating incoming with the outgoing modes, which 
    are indexed in the same way. The first dimension of the array are the outgoing 
    and the second dimension the incoming modes

    Args:
        kzs (float): Z component of the waves
        rhos (float or complex): Mass densities

    Returns:
        complex (2, 2, 1, 1)-array
    """,
    res = np.zeros((2, 2, 1, 1), complex)
    res[1][1][0][0] = 2 * rhos[0] * kzs[1] / (rhos[0] * kzs[1] + rhos[1] * kzs[0])
    res[0][0][0][0] = 2 * rhos[1] * kzs[0] / (rhos[0] * kzs[1] + rhos[1] * kzs[0])
    res[1][0][0][0] = (-rhos[0] * kzs[1] + rhos[1] * kzs[0]) / (rhos[0] * kzs[1] + rhos[1] * kzs[0])
    res[0][1][0][0] = -res[1][0][0][0]
    return res