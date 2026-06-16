"""Scalar spherical wave module.

.. autosummary::
   :toctree:

   periodic_to_scw
   periodic_to_spw
   rotate
   translate
   translate_periodic

"""

import acoustotreams.special as ats
import treams.special as sc
import numpy as np
from treams import lattice
import scipy.special as ss

def _translate_s(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs):
    """Regular translation coefficient for spherical modes"""
    if abs(kr) < 1e-16:
        return 0.0j
    return ats.tl_ssw(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs)

def _translate_r(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs):
    """Singular translation coefficient for spherical modes"""
    return ats.tl_ssw_r(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs)


_translate_s = np.vectorize(_translate_s)
_translate_r = np.vectorize(_translate_r)


def translate(lambda_, mu, l, m, kr, theta, phi, singular=True, *args, **kwargs):
    """translate(lambda_, mu, l, m, kr, theta, phi, singular=True)

    Translation coefficients for spherical waves.

    Returns the correct translation coefficients from :func:acoustotreams.ssw.tl_ssw
    and :func:acoustotreams.ssw.tl_ssw_r or combinations thereof for the specified modes and
    basis.

    Args:
        lambda_ (int, array_like): Degree of output modes.
        mu (int, array_like): Order of output modes.
        l (int, array_like): Degree of input modes.
        m (int, array_like): Order of input modes.
        kr (float or complex, array_like): Translation distance in units of the wavenumber.
        theta (float, array_like): Polar angle (rad).
        phi (float, array_like): Azimuthal angle (rad).
        singular (bool, optional): If True, singular translation coefficients are used,
            else regular coefficients. Defaults to ``True``.

    Returns:
        complex
    """
    if singular:
        return _translate_s(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs)
    return _translate_r(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs)


def _rotate(lambda_, mu, l, m, phi, theta, psi, *args, **kwargs):
    """
    Rotation coefficients for the rotation by the Euler angles phi, theta, and psi.
    It is intended to be used for constructing the rotation matrix manually.
    """
    if lambda_ == l:
        return sc.wignerd(l, mu, m, phi, theta, psi, *args, **kwargs)
    else:
        return 0.0j

_rotate = np.vectorize(_rotate)


def rotate(lambda_, mu, l, m, phi, theta=0, psi=0, *args, **kwargs):
    """rotate(lambda_, mu, l, m, phi, theta=0, psi=0)
    
    Rotation coefficients for scalar spherical modes.

    Return the correct rotation coefficients from :func:`treams.special.wignerd`. The
    angles are given as the Euler angles in `z-y-z`-convention. In the intrinsic (object
    fixed coordinate system) convention, the rotations are applied in the order phi
    first, theta second, and psi third. In the extrinsic (global or reference frame fixed
    coordinate system), the rotations are applied psi first, theta second, phi third.

    Args:
        lambda_ (int, array_like): Degree of output modes.
        mu (int, array_like): Order of output modes.
        l (int, array_like): Degree of input modes.
        m (int, array_like): Order of input modes.
        phi (float or complex, array_like): First Euler angle.
        theta (float, array_like): Second Euler angle.
        psi (float, array_like): Third Euler angle.

    Returns:
        complex
    """
    return _rotate(lambda_, mu, l, m, phi, theta, psi, *args, **kwargs)


def _transl_a_lattice(lambda_, mu, l, m, dlms):
    """Singular translation coefficient for a scalar spherical wave on a lattice"""
    pref = np.power(-1, np.abs(m)) * np.sqrt(4 * np.pi * (2 * l + 1) * (2 * lambda_ + 1)) * np.power(1.0j, lambda_ - l)
    dlm = 0
    res = 0
    for p in range(l + lambda_, max(abs(lambda_ - l), abs(mu - m)) - 1, -2):
        dlm = dlms[p * (p + 1) + m - mu]
        res += dlm * np.power(1.0j, p) * np.sqrt(2 * p + 1) * sc.wigner3j(l, lambda_, p, m, -mu, -m + mu) * sc.wigner3j(l, lambda_, p, 0, 0, 0)
    return res * pref


def translate_periodic(ks, kpar, a, rs, out, in_=None, rsin=None, eta=0, func=lattice.lsumsw):
    """Translation coefficients for scalar spherical waves on a lattice.

    Return the translation coefficents for the given modes on a lattice. The computations
    employ the fast converging Ewald summation from :py:data:`~treams.lattice`.


    Args:
        ks (float or complex): Wavenumber of longitudinal waves in the medium.
        kpar (float, (D,)-array): Parallel component of the wave vector. Defines the dimension with `1 <= D <= 3`.
        a (float, (D,D)-array): Lattice vectors in each row of the array.
        rs (float, (M, 3)-array): Shift vectors with respect to one lattice point.
        out (2- or 3-tuple of integer arrays): Output modes.
        in_ (2- or 3-tuple of integer arrays): Input modes. If None is given, they are equal to
            the output modes.
        rsin (float): Shift vectors of the input modes. If None is given, they are equal to `rs`.
        eta (float or complex, optional): Splitting parameter for the Ewald summation of
            lattice sums. If zero, an estimation for the optimal value is done. Defaults
            to zero.

    Returns:
        complex array
    """
    if in_ is None:
        in_ = out
    out = (*(np.array(o) for o in out),)
    in_ = (*(np.array(i) for i in in_),)
    if len(out) < 2 or len(out) > 3:
        raise ValueError(f"invalid length of output modes {len(out)}, must be 2 or 3")
    elif len(out) == 2:
        out = (np.zeros_like(out[0]),) + out
    if len(in_) < 2 or len(in_) > 3:
        raise ValueError(f"invalid length of input modes {len(in_)}, must be 2 or 3")
    elif len(in_) == 2:
        in_ = (np.zeros_like(in_[0]),) + in_
    if rsin is None:
        rsin = rs
    modes = np.array([
        [l, m]
        for l in range(np.max(out[1]) + np.max(in_[1]) + 1)
        for m in range(-l, l + 1)
    ])
    kpar = np.array(kpar)
    rs = np.array(rs)
    if rs.ndim == 1:
        rs = rs.reshape((1, -1))
    rsin = np.array(rsin)
    if rsin.ndim == 1:
        rsin = rsin.reshape((1, -1))
    rsdiff = -rs[:, None, None, :] + rsin[:, None, :]

    dim = 1 if kpar.ndim == 0 else kpar.shape[-1]
    # The result has the shape (n_rs, n_rs, n_modes)
    dlms = func(dim, modes[:, 0], modes[:, 1], ks, kpar, a, rsdiff, eta)
    res = np.zeros((out[1:][0].shape[-1], in_[1:][0].shape[-1]), complex)
    for i in range(out[1:][0].shape[-1]):
        for j in range(in_[1:][0].shape[-1]):
                res[i][j] = _transl_a_lattice(out[1:][0][i], out[1:][1][i], in_[1:][0][j], in_[1:][1][j], dlms[out[0][i]][in_[0][j]])
    return res
         

def _periodic_to_spw(kx, ky, kz, l, m, area, *args, **kwargs):
    """Convert a periodic spherical wave into plane waves"""
    k = np.sqrt(kx * kx + ky * ky + kz * kz)
    phi = np.arctan2(ky, kx)
    kz_s = kz
    if kz == 0:
        kz_s = 1e-20 + 1e-20j
    elif (np.imag(kz_s) < 0) or (np.imag(kz_s) == 0 and np.real(kz_s) < 0):
        kz_s = -kz_s
    prefactor = (
        np.sqrt(np.pi * (2 * l + 1))
        * np.sqrt(ss.gamma(l - m + 1) / ss.gamma(l + m + 1))
        * np.exp(1j * m * phi)
        * np.power(-1j, l)
        / (area * k * kz_s)
    )
    return prefactor * sc.lpmv(m, l, kz / k + 0j, *args, **kwargs)

_periodic_to_spw = np.vectorize(_periodic_to_spw)


def periodic_to_spw(kx, ky, kz, l, m, area, *args, **kwargs):
    """periodic_to_spw(kx, ky, kz, l, m, area)
    
    Convert periodic scalar spherical waves into plane waves.

    Return the coefficient for the basis change from spherical waves on a lattice
    to plane waves. For multiple positions, only the diagonal values (corresponding to
    identical positions) are returned. A correct phase factor is still necessary for the correct
    result.

    Args:
        kx (float, array_like): X-component of the wave vector of plane waves.
        ky (float, array_like): Y-component of the wave vector of plane waves.
        kz (float or complex, array_like): Z-component of the wave vector of plane waves.
        l (int, array_like): Degree of spherical waves.
        m (int, array_like): Order of spherical waves.
        area (float, array_like): Unit cell area.

    Returns:
        complex
    """
    return _periodic_to_spw(kx, ky, kz, l, m, area, *args, **kwargs)

def _periodic_to_scw(kz, m, l, mu, k, a, *args, **kwargs):
    """Convert periodic spherical wave to cylindrical wave"""
    if mu != m:
        return 0.0j + sc.lpmv(0, 1, 0, *args, **kwargs)
    prefactor = (
        np.sqrt(np.pi * (2 * l + 1), *args, **kwargs) *
        np.sqrt(ss.gamma(l - m + 1) / ss.gamma(l + m + 1), *args, **kwargs) * 
        0.5 * 
        np.power(-1j, l - m, *args, **kwargs)
        / (a * k)
    ) 
    return prefactor * sc.lpmv(m, l, kz/k+0j, *args, **kwargs)

 
_periodic_to_scw = np.vectorize(_periodic_to_scw)


def periodic_to_scw(kz, m, l, mu, k, area, *args, **kwargs):
    """periodic_to_scw(kz, m, l, mu, k, area)
    
    Convert periodic spherical waves into cylindrical waves.

    Return the coefficient for the basis change from spherical waves on a lattice
    to cylindrical waves. For multiple positions, only the diagonal values (corresponding
    to identical positions) are returned. 

    Args:
        kz (float, array_like): Z-component of the wave vector of cylinbdrical waves.
        m (int, array_like): Order or cylinbdrical waves.
        l (int, array_like): Degree of spherical waves.
        mu (int, array_like): Order of spherical waves.
        k (float or complex, array_like): Angular wavenumber.
        area (float, array_like): Unit cell area.

    Returns:
        complex
    """
    return _periodic_to_scw(kz, m, l, mu, k, area, *args, **kwargs)