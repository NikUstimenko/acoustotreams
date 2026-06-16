"""Scalar plane wave module.

.. autosummary::
   :toctree:

   to_scw
   to_ssw
   translate
   permute_xyz

"""

import numpy as np
import treams.special as sc
import scipy.special as ss


def translate(kx, ky, kz, x, y, z, *args, **kwargs):
    r"""translate(kx, ky, kz, x, y, z)
    
    Translation coefficients for scalar plane waves.

    The translation coefficient is the phase factor
    :math:`\mathrm e^{\mathrm i \mathbf k \mathbf r}`.

    Args:
        kx, ky, kz (float or complex, array_like): Wave vector components.
        x, y, z (float, array_like): Translation vector components.

    Returns:
        complex
    """
    return np.exp(1j * (kx * x + ky * y + kz * z), *args, **kwargs)


def _to_ssw(l, m, kx, ky, kz, *args, **kwargs):
    """Coefficient for the expansion of a scalar plane wave in scalar spherical waves"""
    phi = np.arctan2(ky, kx)
    k = np.sqrt(kx * kx + ky * ky + kz * kz)
    pref = (
        2 * np.sqrt(np.pi * (2 * l + 1))
        * np.sqrt(ss.gamma(l - m + 1) / ss.gamma(l + m + 1))
        * np.power(1j, l)
        * np.exp(-1j * m * phi)
    )
    return pref * sc.lpmv(m, l, kz / k, *args, **kwargs)

_to_ssw = np.vectorize(_to_ssw)

def to_ssw(l, m, kx, ky, kz, *args, **kwargs):
    """to_ssw(l, m, kx, ky, kz)
    
    Coefficients for the expansion of scalar plane waves in scalar spherical waves.

    Return the coefficient for the basis change from scalar plane waves to scalar spherical waves.
    For multiple positions, only the diagonal values (corresponding to identical positions) are
    returned.

    Args:
        l (int, array_like): Degree of scalar spherical waves.
        m (int, array_like): Order of scalar spherical waves.
        kx (float, array_like): X-component of the wave vector of scalar plane waves.
        ky (float, array_like): Y-component of the wave vector of scalar plane waves.
        kz (float, array_like): Z-component of the wave vector of scalar plane waves.

    Returns:
        complex
    """
    return _to_ssw(l, m, kx, ky, kz, *args, **kwargs)

    
def _to_scw(kzcw, m, kx, ky, kzpw, *args, **kwargs):
    """Coefficient for the expansion of a scalar plane wave in scalar cylindricrical waves"""  
    krho = np.sqrt(kx * kx + ky * ky)
    if np.abs(kzcw - kzpw) <= 1e-12:
        if m == 0:
            return np.power(1, 0, *args, **kwargs)
        if krho == 0:
            return np.power(1j, m, *args, **kwargs)
        return np.power((1j * kx + ky) / krho, m, *args, **kwargs)
    elif np.abs(kzcw - kzpw) > 1e-12:
        return 0.0j + sc.lpmv(0, 1, 0, *args, **kwargs)
    

_to_scw = np.vectorize(_to_scw)    

def to_scw(kzcw, m, kx, ky, kzpw, *args, **kwargs):
    """to_scw(qz, m, kx, ky, kz)

    Coefficients for the expansion of scalar plane waves in scalar cylindrical waves.

    Return the coefficient for the basis change from scalar plane waves to scalar cylindrical waves.
    For multiple positions, only the diagonal values (corresponding to identical positions) are returned.

    Args:
        qz (float, array_like): Z-component of the wave vector of scalar cylindrical waves.
        m (int, array_like): Order of scalar cylindrical waves.
        kx (float, array_like): X-component of the wave vector of scalar plane waves.
        ky (float, array_like): Y-component of the wave vector of scalar plane waves.
        kz (float, array_like): Z-component of the wave vector of scalar plane waves.

    Returns:
        complex
    """
    return _to_scw(kzcw, m, kx, ky, kzpw, *args, **kwargs)


def _xyz_to_zxy(kx, ky, kz, *args, **kwargs):
    return np.power(1, 0, *args, **kwargs)

_xyz_to_zxy = np.vectorize(_xyz_to_zxy)  

def _xyz_to_yzx(kx, ky, kz, *args, **kwargs):
    return np.power(1, 0, *args, **kwargs)

_xyz_to_yzx = np.vectorize(_xyz_to_yzx)

def permute_xyz(kx, ky, kz, inverse=False, *args, **kwargs):
    """permute_xyz(kx, ky, kz, inverse=False)
    
    Change the coordinate system of the plane-wave basis.

    A plane wave in the coordinate system :math:`(x, y, z)` with the primary direction of
    propagation along the z-axis is described in the system :math:`(x', y', z') = (y, z, x)`. 
    The inverse transformation is also possible.

    The function is essentially diagonal relative to the wavenumber, because we always describe
    the input and output modes in the unprimed coordinate system.

    Args:
        kxa (float, array_like): X-component of the wave vector of output modes.
        kya (float or complex, array_like): Y-component of the wave vector of output modes.
        kza (float, array_like): Z-component of the wave vector of output modes.
        kx (float, array_like): X-component of the wave vector of input modes.
        ky (float, array_like): Y-component of the wave vector of input modes.
        kz (float or complex, array_like): Z-component of the wave vector of input modes.
        inverse (bool, optional): Defaults to False. If True, use the inverse transformation.

    Returns:
        complex
    """
    if inverse:
        return _xyz_to_zxy(kx, ky, kz, *args, **kwargs)
    return _xyz_to_yzx(kx, ky, kz, *args, **kwargs)