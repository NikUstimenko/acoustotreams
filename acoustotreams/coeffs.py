"""Scattering coefficients for high-symmetry cases.

Calculate the scattering coefficients for cases where they can be obtained analytically
easily. This is a (multilayered) sphere using spherical waves (Mie coefficients), 
a cylinder using cylindrical waves, and an infinitely extended planar
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

def _mie_acoustics_iter(tm, l, x, mat_sphere, mat_env):
    """Solve a matrix for the boundary conditions at a spherical interface.

    Six types of interfaces are supported: fluid-fluid, fluid-solid, solid-fluid,
    solid-solid, soft-fluid, and hard-fluid.

    The fluid medium has :math:`(\rho, c, c_t = 0)`, the solid medium has
    :math:`(\rho, c, c_t)`, the soft medium has :math:`(\rho = 0, c = 0, c_t = 0)`,
    and the hard medium has :math:`(\rho = \infty, c = 0, c_t = 0)`.

    Note:
        The order of coefficients is LL, NL, LN, and NN.

    Args:
        tm: T-matrix of the previous layer.
        l: Degree of the coefficient.
        x: Size parameters in the air.
        mat_sphere: Material of the sphere :math:`(\rho, c, c_t)`.
        mat_env: Material of the medium :math:`(\rho, c, c_t)`.

    Returns:
        complex 4-array
    """ 
    x_env = x * AcousticMaterial().c / mat_env[1]

    j = treams.special.spherical_jn(l, x_env, derivative=False)
    h = treams.special.spherical_hankel1(l, x_env)
    if (np.abs(mat_sphere[2]) == 0 
        and np.abs(mat_sphere[1]) == 0 
        and np.abs(mat_sphere[0]) == 0
        and np.abs(mat_env[2]) == 0):
        return np.array([-j / h, 0, 0, 0])

    j_d = treams.special.spherical_jn(l, x_env, derivative=True)
    h_d = treams.special.spherical_hankel1_d(l, x_env)
    if (np.abs(mat_sphere[2]) == 0 
        and np.abs(mat_sphere[1]) == 0 
        and np.abs(mat_sphere[0]) == np.inf
        and np.abs(mat_env[2]) == 0):
        return np.array([-j_d / h_d, 0, 0, 0])
    
    x_sphere = x * AcousticMaterial().c / mat_sphere[1]
    j1 = treams.special.spherical_jn(l, x_sphere, derivative=False)
    j1_d = treams.special.spherical_jn(l, x_sphere, derivative=True)
    h1 = treams.special.spherical_hankel1(l, x_sphere)
    h1_d = treams.special.spherical_hankel1_d(l, x_sphere)

    psi = np.sqrt(l * (l + 1))
    d12 = h
    d22 = x_env * h_d
    d32 = x_env * h_d - h
    d14 = j1
    f14 = h1
    d24 = x_sphere * j1_d
    f24 = x_sphere * h1_d
    d1L = j
    d2L = x_env * j_d
    d3L = x_env * j_d - j
    
    res = np.zeros(4, complex)

    if np.abs(mat_env[2]) > 0:
        x_env_t = x * AcousticMaterial().c / mat_env[2]
        d42 = (psi**2 - 0.5 * x_env_t**2) * h - 2 * x_env * h_d
        d4L = (psi**2 - 0.5 * x_env_t**2) * j - 2 * x_env * j_d
        if l > 0:
            jt = treams.special.spherical_jn(l, x_env_t, derivative=False)
            jt_d = treams.special.spherical_jn(l, x_env_t, derivative=True)
            ht = treams.special.spherical_hankel1(l, x_env_t)
            ht_d = treams.special.spherical_hankel1_d(l, x_env_t)
            d21 = psi**2 * ht 
            d31 = (psi**2 - 0.5 * x_env_t**2 - 1) * ht - x_env_t * ht_d
            d41 = psi**2 * (x_env_t * ht_d - ht)
            d2N = psi**2 * jt
            d3N = (psi**2 - 0.5 * x_env_t**2 - 1) * jt - x_env_t * jt_d
            d4N = psi**2 * (x_env_t * jt_d - jt) 

    if np.abs(mat_sphere[2]) > 0:
        x_sphere_t = x * AcousticMaterial().c / mat_sphere[2]
        if l > 0:
            j1t = treams.special.spherical_jn(l, x_sphere_t, derivative=False)
            j1t_d = treams.special.spherical_jn(l, x_sphere_t, derivative=True)
            h1t = treams.special.spherical_hankel1(l, x_sphere_t)
            h1t_d = treams.special.spherical_hankel1_d(l, x_sphere_t)
            f23 = psi**2 * h1t
            d23 = psi**2 * j1t

    if np.abs(mat_sphere[2]) == 0 and np.abs(mat_env[2]) == 0:
        delta = x_sphere * mat_env[0] / (x_env * mat_sphere[0]) 
        matrix = np.zeros((2, 2), complex)
        rhs = np.zeros(2, complex)
        matrix[0][0] = -d12
        matrix[0][1] = (d14 + tm[0] * f14) / delta
        rhs[0] = d1L
        matrix[1][0] = -d22 / x_env
        matrix[1][1] = (d24 + tm[0] * f24) / x_sphere
        rhs[1] = d2L / x_env
        sol = np.linalg.solve(matrix, rhs)
        res[0] = sol[0]
        return res
    elif np.abs(mat_sphere[2]) == 0 and np.abs(mat_env[2]) > 0:
        d44 = -mat_sphere[0] / mat_env[0] * 0.5 * x_env_t**2 * j1 
        f44 = -mat_sphere[0] / mat_env[0] * 0.5 * x_env_t**2 * h1 
        if l == 0:
            matrix = np.zeros((2, 2), complex)
            rhs = np.zeros(2, complex)
            matrix[0][0] = -d22 / x_env
            matrix[0][1] = (d24 + tm[0] * f24) / x_sphere
            matrix[1][0] = -d42 / x_env
            matrix[1][1] = (d44 + tm[0] * f44) / x_sphere
            rhs[0] = d2L / x_env
            rhs[1] = d4L / x_env
            sol = np.linalg.solve(matrix, rhs)
            res[0] = sol[0]
            return res
        elif l > 0:
            matrix = np.zeros((3, 3), complex)
            rhs_L = np.zeros(3, complex)
            rhs_N = np.zeros(3, complex)
            matrix[0][0] = d21 / x_env_t
            matrix[0][1] = -psi * d22 / x_env
            matrix[0][2] = psi * (d24 + tm[0] * f24) / x_sphere
            rhs_L[0] = psi * d2L / x_env
            rhs_L[1] = psi * d3L / x_env
            rhs_L[2] = psi * d4L / x_env
            matrix[1][0] = d31 / x_env_t
            matrix[1][1] = -psi * d32 / x_env
            matrix[2][0] = d41 / x_env_t
            matrix[2][1] = -psi * d42 / x_env
            matrix[2][2] = psi * (d44 + tm[0] * f44) / x_sphere
            sol_L = np.linalg.solve(matrix, rhs_L)
            res[0] = sol_L[1]
            res[1] = sol_L[0]
            rhs_N[0] = -d2N / x_env
            rhs_N[1] = -d3N / x_env
            rhs_N[2] = -d4N / x_env
            sol_N = np.linalg.solve(matrix, rhs_N)
            res[2] = sol_N[1]
            res[3] = sol_N[0]
            return res
    elif np.abs(mat_sphere[2]) > 0 and np.abs(mat_env[2]) == 0:
        d42 = -mat_env[0]/mat_sphere[0] * x_sphere_t**2 * 0.5 * h
        d44 = (l * (l + 1) - 0.5 * x_sphere_t**2) * j1 - 2 * x_sphere * j1_d
        f44 = (l * (l + 1) - 0.5 * x_sphere_t**2) * h1 - 2 * x_sphere * h1_d
        d4L = -mat_env[0]/mat_sphere[0] * 0.5 * x_sphere_t**2 * j
        if l == 0:
            matrix = np.zeros((2, 2), complex)
            rhs = np.zeros(2, complex)
            matrix[0][0] = -d22 / x_env
            matrix[0][1] = (d24 + tm[0] * f24) / x_sphere
            matrix[1][0] = -d42 / x_env
            matrix[1][1] = (d44 + tm[0] * f44) / x_sphere
            rhs[0] = d2L / x_env
            rhs[1] = d4L / x_env
            sol = np.linalg.solve(matrix, rhs)
            res[0] = sol[0]
            return res
        elif l > 0:
            d33 = (l * (l + 1) - 0.5 * x_sphere_t**2 - 1) * j1t - x_sphere_t * j1t_d
            d43 = l * (l + 1) * (x_sphere_t * j1t_d - j1t)
            f33 = (psi**2 - 0.5 * x_sphere_t**2 - 1) * h1t - x_sphere_t * h1t_d
            f43 = psi**2 * (x_sphere_t * h1t_d - h1t)  
            d34 = x_sphere * j1_d - j1    
            f34 = x_sphere * h1_d - h1 
            matrix = np.zeros((3, 3), complex)
            rhs = np.zeros(3, complex)
            matrix[0][0] = -psi * d22 / x_env
            matrix[2][0] = -psi * d42 / x_env
            matrix[0][1] = -d23 / x_sphere_t - f23 * tm[3] / x_sphere_t + psi * tm[2] * f24 / x_sphere
            matrix[1][1] = -d33 / x_sphere_t - f33 * tm[3] / x_sphere_t + psi * tm[2] * f24 / x_sphere
            matrix[2][1] = -d43 / x_sphere_t - f43 * tm[3] / x_sphere_t + psi * tm[2] * f24 / x_sphere
            matrix[0][2] = psi * (d24 + tm[0] * f24) / x_sphere - psi * tm[1] * f23 / x_sphere_t
            matrix[1][2] = psi * (d34 + tm[0] * f34) / x_sphere - psi * tm[1] * f33 / x_sphere_t
            matrix[2][2] = psi * (d44 + tm[0] * f44) / x_sphere - psi * tm[1] * f43 / x_sphere_t 
            rhs[0] = psi * d2L / x_env
            rhs[2] = psi * d4L / x_env
            sol = np.linalg.solve(matrix, rhs)
            res[0] = sol[0]
            return res 
    elif np.abs(mat_sphere[2]) > 0 and np.abs(mat_env[2]) > 0:
        d44 = (mat_sphere[0]/mat_env[0] * (x_env_t/x_sphere_t)**2 
               * ((psi**2 - 0.5 * x_sphere_t**2) * j1 - 2 * x_sphere * j1_d))
        f44 = (mat_sphere[0]/mat_env[0] * (x_env_t/x_sphere_t)**2 
               * ((psi**2 - 0.5 * x_sphere_t**2) * h1 - 2 * x_sphere * h1_d))
        if l == 0:
            matrix = np.zeros((2, 2), complex)
            rhs_L = np.zeros(2, complex)
            matrix[0][0] = -d22 / x_env
            matrix[1][0] = -d42 / x_env
            matrix[0][1] = (d24 + tm[0] * f24) / x_sphere
            matrix[1][1] = (d44 + tm[0] * f44) / x_sphere
            rhs_L[0] = d2L / x_env
            rhs_L[1] = d4L / x_env
            sol_L = np.linalg.solve(matrix, rhs_L)
            res[0] = sol_L[0]
            return res
        elif l > 0:
            d11 = x_env_t * ht_d + ht
            d13 = x_sphere_t * j1t_d + j1t
            d33 = (mat_sphere[0]/mat_env[0] * (x_env_t/x_sphere_t)**2 
                * ((psi**2 - 0.5 * x_sphere_t**2 - 1) * j1t - x_sphere_t * j1t_d))
            d43 = (mat_sphere[0]/mat_env[0] * (x_env_t/x_sphere_t)**2  * psi**2 
                * (x_sphere_t * j1t_d - j1t))
            d34 = (mat_sphere[0]/mat_env[0] * (x_env_t/x_sphere_t)**2 
                * (x_sphere * j1_d - j1))
            f34 = (mat_sphere[0]/mat_env[0] * (x_env_t/x_sphere_t)**2 
               * (x_sphere * h1_d - h1))
            f13 = x_sphere_t * h1_d + h1
            f33 = (mat_sphere[0]/mat_env[0] * (x_env_t/x_sphere_t)**2 
               * ((psi**2 - 0.5 * x_sphere_t**2 - 1) * h1t - x_sphere_t * h1t_d))
            f43 = (mat_sphere[0]/mat_env[0] * (x_env_t/x_sphere_t)**2  * psi**2 
               * (x_sphere_t * h1t_d - h1t))
            d1N = x_env_t * jt_d + jt
            matrix = np.zeros((4, 4), complex)
            rhs_L = np.zeros(4, complex)
            rhs_N = np.zeros(4, complex)
            matrix[0][0] = d11 / x_env_t
            matrix[1][0] = d21 / x_env_t
            matrix[2][0] = d31 / x_env_t
            matrix[3][0] = d41 / x_env_t
            matrix[0][1] = -psi * d12 / x_env
            matrix[1][1] = -psi * d22 / x_env
            matrix[2][1] = -psi * d32 / x_env
            matrix[3][1] = -psi * d42 / x_env
            matrix[0][2] = -d13 / x_sphere_t - tm[3] * f13 / x_sphere_t + psi * tm[2] * f14 / x_sphere
            matrix[1][2] = -d23 / x_sphere_t - tm[3] * f23 / x_sphere_t + psi * tm[2] * f24 / x_sphere
            matrix[2][2] = -d33 / x_sphere_t - tm[3] * f33 / x_sphere_t + psi * tm[2] * f34 / x_sphere
            matrix[3][2] = -d43 / x_sphere_t - tm[3] * f43 / x_sphere_t + psi * tm[2] * f44 / x_sphere
            matrix[0][3] = psi * (d14 + tm[0] * f14) / x_sphere - tm[1] * f13 / x_sphere_t 
            matrix[1][3] = psi * (d24 + tm[0] * f24) / x_sphere - tm[1] * f23 / x_sphere_t 
            matrix[2][3] = psi * (d34 + tm[0] * f34) / x_sphere - tm[1] * f33 / x_sphere_t 
            matrix[3][3] = psi * (d44 + tm[0] * f44) / x_sphere - tm[1] * f43 / x_sphere_t 
            rhs_L[0] = psi * d1L / x_env
            rhs_L[1] = psi * d2L / x_env
            rhs_L[2] = psi * d3L / x_env
            rhs_L[3] = psi * d4L / x_env
            rhs_N[0] = -d1N / x_env_t
            rhs_N[1] = -d2N / x_env_t
            rhs_N[2] = -d3N / x_env_t
            rhs_N[3] = -d4N / x_env_t
            sol_N = np.linalg.solve(matrix, rhs_N)
            res[2] = sol_N[1]
            res[3] = sol_N[0] 
            sol_L = np.linalg.solve(matrix, rhs_L)
            res[0] = sol_L[1]
            res[1] = sol_L[0] 
            return res
    return res 

def mie_acoustics(l, x, *materials):
    r"""Mie scattering coefficient of degree l.

    The sphere is defined by its size parameter :math:`k_0 r`, where :math:`r` is the
    radius and :math:`k_0` the wave number in air. A multilayered sphere is defined
    by giving an array of ascending numbers, that define the size parameters of the
    sphere and its shells starting from the center.

    Likewise, the material parameters are given from inside to outside. These arrays
    are expected to be exactly one unit larger then the array `x`.

    The result is a complex number relating incident with the scattered modes, which are 
    index in the same way.

    Note:
        1. :math:`c_t` of the last material must be zero. 
        2. For the soft and hard spheres, only one radius must be given. 

    Args:
        l (integer): Degree :math:`l \geq 0`
        x (float, array_like): Size parameters
        rho (float or complex, array_like): Mass density
        c (float or complex, array_like): Longitudinal speed of sound
        c_t (float or complex, array_like): Transverse speed of sound

    Returns:
        complex array
    """
    
    mat = list(zip(*materials))
    l = np.atleast_1d(l)
    x = np.atleast_1d(x)
    res = []
    for j in l:
        mie = np.zeros(4, complex)
        for i in range(len(list(mat)) - 1):
            mat_sphere, mat_env = mat[i], mat[i + 1]
            mie = _mie_acoustics_iter(mie, j, x[i], mat_sphere, mat_env)
        res.append(mie[0])
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


#if (np.abs(mat_sphere[2]) == 0 
#    and np.abs(mat_sphere[1]) == 0 
#    and np.abs(mat_sphere[0]) == 0
#    and np.abs(mat_env[2]) > 0):
#    if l == 0:
#        return np.array([-j / h, 0, 0, 0])
#    matrix = np.zeros((2, 2), complex)
#    rhs_L = np.zeros(2, complex)
#    rhs_N = np.zeros(2, complex)
#    matrix[0][0] = d31 / x_env_t
#    matrix[0][1] = -psi * d32 / x_env
#    matrix[1][0] = d41 / x_env_t
#    matrix[1][1] = -psi * d42 / x_env
#    rhs_L[0] = psi * d3L / x_env
#    rhs_L[1] = psi * d4L / x_env
#    rhs_N[0] = -d3N / x_env_t
#    rhs_N[1] = -d4N / x_env_t
#    sol_L = np.linalg.solve(matrix, rhs_L)
#    res[0] = sol_L[1]
#    res[1] = sol_L[0]
#    sol_N = np.linalg.solve(matrix, rhs_N)
#    res[2] = sol_N[1]
#    res[3] = sol_N[0]
#    return res

#if (np.abs(mat_sphere[2]) == 0 
#    and np.abs(mat_sphere[1]) == 0 
#    and np.abs(mat_sphere[0]) == np.inf
#    and np.abs(mat_env[2]) > 0):
#    if l == 0:
#        return np.array([-j_d / h_d, 0, 0, 0])
#    matrix = np.zeros((2, 2), complex)
#    rhs_L = np.zeros(2, complex)
#    rhs_N = np.zeros(2, complex)
#    matrix[0][0] = d11 / x_env_t
#    matrix[0][1] = -psi * d12 / x_env
#    matrix[1][0] = d21 / x_env_t
#    matrix[1][1] = -psi * d22 / x_env
#    rhs_L[0] = psi * d1L / x_env
#    rhs_L[1] = psi * d2L / x_env
#    rhs_N[0] = -d1N / x_env_t
#    rhs_N[1] = -d2N / x_env_t
#    sol_L = np.linalg.solve(matrix, rhs_L)
#    res[0] = sol_L[1]
#    res[1] = sol_L[0]
#    sol_N = np.linalg.solve(matrix, rhs_N)
#    res[2] = sol_N[1]
#    res[3] = sol_N[0]
#    return res