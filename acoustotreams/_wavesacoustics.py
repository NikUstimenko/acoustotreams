import treams.special as sc
import numpy as np
import cmath
import scipy.special as ss

def spw_Psi(kx, ky, kz, x, y, z):
    """Scalar plane wave"""
    return np.exp(1j * (kx * x + ky * y + kz * z))


def vpw_L(kx, ky, kz, x, y, z):
    """Longitudinal vector plane wave"""
    k = np.sqrt(kx * kx + ky * ky + kz * kz)
    phase = np.exp(1j * (kx * x + ky * y + kz * z))
    if k == 0:
        return np.asarray([np.nan, np.nan, np.nan]).T
    else:
        return 1j * np.asarray([kx * phase, ky * phase, kz * phase], np.complex128).T / k

def ssw_Psi(l, m, kr, theta, phi):
     """Singular scalar spherical wave"""
     return sc.sph_harm(m, l, phi, theta) * sc.spherical_hankel1(l, kr)

def ssw_psi(l, m, ksx, ksy, ksz, theta, phi, ks):
     """Far-field amplitude of singular scalar spherical wave"""
     return (sc.sph_harm(m, l, phi, theta)
             * np.exp(-1j * ((ksx * np.cos(phi) + ksy * np.sin(phi)) * np.sin(theta) 
                             + ksz * np.cos(theta)))
             * np.power(-1j, l + 1)
             * 1 / ks)

def ssw_rPsi(l, m, kr, theta, phi):
     """Regular scalar spherical wave"""
     return sc.sph_harm(m, l, phi, theta) * sc.spherical_jn(l, kr)


def scw_Psi(kz, m, krr, phi, z):
     """Singular scalar cylindrical wave"""
     return np.exp(1j * (m * phi + kz * z)) * sc.hankel1(m, krr)

def scw_psi(kz, m, krhox, krhoy, phi, z, krho):
     """Far-field amplitude of singular scalar cylindrical wave"""
     return (np.exp(1j * m * phi + 1j * kz * z) * 
             np.sqrt(2 / (np.pi * krho)) *
             np.exp(-1j * ((krhox * np.cos(phi) + krhoy * np.sin(phi)))) *
             np.power(-1j, m) * 
             np.exp(-1j * np.pi/4))

def scw_rPsi(kz, m, krr, phi, z):
     """Regular scalar cylindrical wave"""
     return np.exp(1j * (m * phi + kz * z)) * sc.jv(m, krr)


def vsw_L(l, m, kr, theta, phi):
     """Longitudinal singular vector spherical wave"""
     return np.transpose(sc.vsh_Z(l, m, theta, phi).T * sc.spherical_hankel1_d(l, kr)  + np.sqrt(
         l * (l + 1)
     ) * sc.vsh_Y(l, m, theta, phi).T * sc.spherical_hankel1(l, kr) / kr) * (-1.0j)


def vsw_l(l, m, ksx, ksy, ksz, theta, phi, ks):
     """Far-field amplitude of longitudinal singular vector spherical wave"""
     return np.transpose(
         np.array(
             [
                 sc.sph_harm(m, l, phi, theta)
                 * np.exp(-1j * ((ksx * np.cos(phi) + ksy * np.sin(phi)) * np.sin(theta) 
                                 + ksz * np.cos(theta)))
                 * np.power(-1j, l)
                 * 1 / ks, 
                 np.zeros(len(l)).T,
                 np.zeros(len(l)).T
            ]
        )
    )


def vsw_rL(l, m, kr, theta, phi):
     """Longitudinal regular vector spherical wave"""
     res = []
     for i, x in enumerate(kr):
            if x == 0:
                if l[i] != 1:
                    val = np.zeros(3, complex).T
                    res.append(val)
                else:
                   val = np.transpose(
                       sc.vsh_Z(1, m[i], 0, 0).T + 
                       np.sqrt(2) * sc.vsh_Y(1, m[i], 0, 0).T) * (1/3) * (-1.0j)
                   res.append(val)
            else:      
                val = np.transpose(
                    sc.vsh_Z(l[i], m[i], theta[i], phi[i]).T * 
                    sc.spherical_jn_d(l[i], kr[i])  
                    + np.sqrt(l[i] * (l[i] + 1)) * 
                    sc.vsh_Y(l[i], m[i], theta[i], phi[i]).T * 
                    sc.spherical_jn(l[i], kr[i]) / kr[i]) * (-1.0j)
                res.append(val)
     return np.array(res)


def vcw_L(kz, m, krr, phi, z, krho, k):
     """Longitudinal singular vector cylindrical wave"""
     return np.transpose(
         np.array(
             [
                 sc.hankel1_d(m, krr) * krho / k  * np.exp(1j * (m * phi + kz * z)),
                 1j * m * sc.hankel1(m, krr) / krr * krho / k  * np.exp(1j * (m * phi + kz * z)),
                 1j * kz / k * sc.hankel1(m, krr) * np.exp(1j * (m * phi + kz * z)),
             ]
         )
     )

def vcw_l(kz, m, krhox, krhoy, phi, z, krho, k):
     """Far-field amplitude of longitudinal singular vector cylindrical wave"""
     return np.transpose(
         np.array(
             np.exp(-1j * np.pi/4) *
             [
                np.sqrt(2 / (np.pi * krho)) *
                np.exp(-1j * ((krhox * np.cos(phi) + krhoy * np.sin(phi)))) *
                np.power(-1j, m) * 
                np.exp(1j * (m * phi + kz * z)) * 
                1j * krho / k,
                np.zeros(len(krho)).T,
                np.sqrt(2 / (np.pi * krho)) *
                np.exp(-1j * ((krhox * np.cos(phi) + krhoy * np.sin(phi)))) *
                np.power(-1j, m) * 
                np.exp(1j * (m * phi + kz * z)) * 
                1j * kz / k,
             ]
         )
     )

def vcw_rL(kz, m, krr, phi, z, krho, k):
     """Longitudinal regular vector cylindrical wave"""
     res = []
     for i, x in enumerate(krr):
            if x == 0:
                if (abs(m[i]) != 1) and (m[i] != 0):
                    val = np.zeros(3, complex).T
                    res.append(val)
                elif abs(m[i]) == 1:
                   val = np.array(
                                  [
                       0.5 * krho[i] / k[i],
                       0.5 * 1j * m[i] * krho[i] / k[i],
                       0. 
                       ]).T * np.exp(1j * (m[i] * phi[i] + kz[i] * z[i]))
                   res.append(val)
                elif m[i] == 0:
                    val = np.array(
                        [
                            0., 
                            0., 
                            1j * kz[i] / k[i] * np.exp(1j * (m[i] * phi[i] + kz[i] * z[i]))
                            ]
                            ).T
                    res.append(val)
            else:      
                val = np.array(
                [
                    sc.jv_d(m[i], krr[i]) * krho[i] / k[i],
                    1j * m[i] * sc.jv(m[i], krr[i]) / krr[i] * krho[i] / k[i],
                    1j * kz[i] / k[i] * sc.jv(m[i], krr[i]),
                ] 
            ).T * np.exp(1j * (m[i] * phi[i] + kz[i] * z[i]))
                res.append(val)
     return np.array(res)


def _tl_ssw_helper(l, m, lambda_, mu, p, q):
    """Helper function for the translation coefficient of scalar and longitudinal spherical waves"""
    if (
        p < max(abs(m + mu), abs(l - lambda_))
        or p > abs(l + lambda_)
        or q < abs(l - lambda_)
        or q > abs(l + lambda_)
        or (q + l + lambda_) % 2 != 0
    ):
        return 0
    return (
        (2 * p + 1)
        * np.power(1j, lambda_ - l + p)
        * np.sqrt(ss.gamma(p - m - mu + 1) / ss.gamma(p + m + mu + 1))
        * sc.wigner3j(l, lambda_, p, m, mu, -(m + mu))
        * sc.wigner3j(l, lambda_, q, 0, 0, 0)
    )


def tl_ssw(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs):
    """Singular translation coefficient of scalar and longitudinal spherical waves"""
    pref = np.power(-1, np.abs(m)) * np.sqrt((2 * l + 1) * (2 * lambda_ + 1)) * cmath.exp(1j * (m - mu) * phi)
    res = 0.
    max_ = np.max([np.abs(int(lambda_) - int(l)), np.abs(int(m) - int(mu))])
    min_ = int(l) +  int(lambda_)
    for p in range(min_, max_ - 1, -2):
        res += (
            _tl_ssw_helper(l, m, lambda_, -mu, p, p)
            * sc.spherical_hankel1(p, kr)
            * sc.lpmv(m - mu, p, np.cos(theta), *args, **kwargs)
        )
    return res * pref


def tl_ssw_r(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs):
    """Regular translation coefficient of scalar and longitudinal spherical waves"""
    pref = np.power(-1, np.abs(m)) * np.sqrt((2 * l + 1) * (2 * lambda_ + 1)) * cmath.exp(1j * (m - mu) * phi)
    res = 0.
    max_ = np.max([np.abs(int(lambda_) - int(l)), np.abs(int(m) - int(mu))])
    min_ = int(l) + int(lambda_)
    for p in range(min_, max_ - 1, -2):
        res += (
            _tl_ssw_helper(l, m, lambda_, -mu, p, p)
            * sc.spherical_jn(p, kr)
            * sc.lpmv(m - mu, p, np.cos(theta), *args, **kwargs)
        )
    return res * pref

def tl_scw(kz1, mu, kz2, m, krr, phi, z, *args, **kwargs):
    """Singular translation coefficient of scalar and longitudinal cylindrical waves"""
    if kz1 != kz2:
        return 0.+0.j
    return sc.hankel1(m - mu, krr) * np.exp(1j * ((m - mu) * phi + kz1 * z), *args, **kwargs)


def tl_scw_r(kz1, mu, kz2, m, krr, phi, z, *args, **kwargs):
    """Regular translation coefficient of scalar and longitudinal cylindrical waves"""
    if kz1 != kz2:
        return 0.+0.j
    return sc.jv(m - mu, krr) * np.exp(1j * ((m - mu) * phi + kz1 * z), *args, **kwargs)