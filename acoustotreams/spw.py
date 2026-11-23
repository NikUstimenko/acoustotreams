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


def translate(kx, ky, kz, x, y, z, where=True):
    r"""translate(kx, ky, kz, x, y, z)
    
    Translation coefficients for scalar plane waves.

    The translation coefficient is the phase factor
    :math:`\mathrm e^{\mathrm i \mathbf k \mathbf r}`.

    Args:
        kx (float or complex or array_like): X-component of the wave vector.
        ky (float or complex or array_like): Y-component of the wave vector.
        kz (float or complex or array_like): Z-component of the wave vector.
        x (float or array_like): X-coordinate. Has units of 1/kx.
        y (float or array_like): Y-coordinate. Has units of 1/ky.
        z (float or array_like): Z-coordinate. Has units of 1/kz.

    Returns:
        complex or array_like
    """
    phase = 1j * (kx * x + ky * y + kz * z)
    if where is True:
        return np.exp(phase)
    # Masked entries are set to zero instead of being left uninitialized
    out = np.zeros(np.broadcast_shapes(np.shape(phase), np.shape(where)), complex)
    return np.exp(phase, out=out, where=where)[()]


def _to_ssw(l, m, kx, ky, kz, where=True):
    """Coefficient for the expansion of a scalar plane wave in scalar spherical waves"""
    krho = np.sqrt(kx * kx + ky * ky)
    # exp(i phi) of a possibly complex wave vector (lossy media); for krho = 0, only
    # m = 0 contributes, so any value can be used
    eiphi = 1 + 0j if krho == 0 else (kx + 1j * ky) / krho
    k = np.sqrt(krho * krho + kz * kz)
    pref = (
        2 * np.sqrt(np.pi * (2 * l + 1))
        * np.sqrt(ss.gamma(l - m + 1) / ss.gamma(l + m + 1))
        * np.power(1j, l)
        * np.power(eiphi, -m)
    )
    return pref * sc.lpmv(m, l, kz / k)

_to_ssw = np.vectorize(_to_ssw, otypes=[complex])

def to_ssw(l, m, kx, ky, kz, *args, **kwargs):
    """to_ssw(l, m, kx, ky, kz)
    
    Coefficients for the expansion of scalar plane waves in scalar spherical waves.

    Return the coefficient for the basis change from scalar plane waves to scalar spherical waves.
    For multiple positions, only the diagonal values (corresponding to identical positions) are
    returned.

    Args:
        l (int or array_like): Degree of scalar spherical waves.
        m (int or array_like): Order of scalar spherical waves.
        kx (float or array_like): X-component of the wave vector of scalar plane waves.
        ky (float or array_like): Y-component of the wave vector of scalar plane waves. Has units of kx.
        kz (float or array_like): Z-component of the wave vector of scalar plane waves. has units of kx.
 
    Returns:
        complex or array_like
    """
    return _to_ssw(l, m, kx, ky, kz, *args, **kwargs)

    
def _to_scw(kzcw, m, kx, ky, kzpw, where=True):
    """Coefficient for the expansion of a scalar plane wave in scalar cylindrical waves"""
    if not where or np.abs(kzcw - kzpw) > 1e-12:
        return 0j
    if m == 0:
        return 1 + 0j  
    krho = np.sqrt(kx * kx + ky * ky)
    if krho == 0:
        return np.power(1j, m)
    return np.power((1j * kx + ky) / krho, m)
    

_to_scw = np.vectorize(_to_scw, otypes=[complex])    

def to_scw(kzcw, m, kx, ky, kzpw, *args, **kwargs):
    """to_scw(qz, m, kx, ky, kz)

    Coefficients for the expansion of scalar plane waves in scalar cylindrical waves.

    Return the coefficient for the basis change from scalar plane waves to scalar cylindrical waves.
    For multiple positions, only the diagonal values (corresponding to identical positions) are returned.

    Args:
        kzcw (float or array_like): Z-component of the wave vector of scalar cylindrical waves.
        m (int or array_like): Order of scalar cylindrical waves.
        kx (float or array_like): X-component of the wave vector of scalar plane waves. Has units of qz.
        ky (float or array_like): Y-component of the wave vector of scalar plane waves. Has units of qz.
        kzpw (float or array_like): Z-component of the wave vector of scalar plane waves. Has units of qz.

    Returns:
        complex or array_like
    """
    return _to_scw(kzcw, m, kx, ky, kzpw, *args, **kwargs)


def _xyz_to_zxy(kx, ky, kz, where=True):
    return 1.0 if where else 0.0

_xyz_to_zxy = np.vectorize(_xyz_to_zxy)  

def _xyz_to_yzx(kx, ky, kz, where=True):
    return 1.0 if where else 0.0

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
        kx (float or array_like): X-component of the wave vector.
        ky (float or array_like): Y-component of the wave vector.
        kz (float or complex or array_like): Z-component of the wave vector.
        inverse (bool, optional): Defaults to False. If True, use the inverse transformation.

    Returns:
        complex or array_like
    """
    if inverse:
        return _xyz_to_zxy(kx, ky, kz, *args, **kwargs)
    return _xyz_to_yzx(kx, ky, kz, *args, **kwargs)