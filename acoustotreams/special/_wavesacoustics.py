"""Special mathematical functions related to the scalar waves"""

import treams.special as sc
import numpy as np
import cmath

def spw_Psi(kx, ky, kz, x, y, z):
    r"""Scalar plane wave :math:`\Psi` 

    The scalar plane wave is defined as

    .. math::

        \Psi_{\mathbf k}(\mathbf r)
        = \mathrm e^{\mathrm i \mathbf k \mathbf r}.

    This function describes a solution to the scalar Helmholtz wave
    equation in Cartesian coordinates.

    Args:
        kx (float or complex or array_like): X-component of the wave vector.
        ky (float or complex or array_like): Y-component of the wave vector.
        kz (float or complex or array_like): Z-component of the wave vector.
        x (float or array_like): X-coordinate. Has units of 1/kx.
        y (float or array_like): Y-coordinate. Has units of 1/ky.
        z (float or array_like): Z-coordinate. Has units of 1/kz.

    Returns:
        complex
    
    """
    return np.exp(1j * (kx * x + ky * y + kz * z))


def vpw_L(kx, ky, kz, x, y, z):
    r"""Longitudinal vector plane wave L
    
    The longitudinal vector plane wave is defined as

    .. math::

       \mathbf L_{\mathbf k}(\mathbf r)
        = \mathrm i \frac{\mathbf k}{k} \mathrm e^{\mathrm i \mathbf k \mathbf r}.

    This function describes a longitudinal solution to the vector Helmholtz wave
    equation in Cartesian coordinates.

    Args:
        kx (float or complex or array_like): X-component of the wave vector.
        ky (float or complex or array_like): Y-component of the wave vector. Has units of kx.
        kz (float or complex or array_like): Z-component of the wave vector. Has units of kx.
        x (float or array_like): X-coordinate. Has units of 1/kx.
        y (float or array_like): Y-coordinate. Has units of 1/kx.
        z (float or array_like): Z-coordinate. Has units of 1/kx.

    Returns:
        complex array
    
    """
    k = np.sqrt(np.asarray(kx * kx + ky * ky + kz * kz))
    phase = np.exp(1j * (kx * x + ky * y + kz * z))
    vec = np.stack(np.broadcast_arrays(kx * phase, ky * phase, kz * phase), axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        res = 1j * vec / k[..., None]
    return np.where((k == 0)[..., None], np.nan, res)

def ssw_Psi(l, m, kr, theta, phi):
    r"""Singular scalar spherical wave :math:`\Psi`

    The singular scalar spherical wave is defined as

    .. math::

        \Psi_{lm}^{(3)}(kr, \theta, \varphi)
        = h_l^{(1)}(kr) Y_{lm}(\theta, \varphi)

    with :py:obj:`treams.special.sph_harm` and :py:obj:`treams.special.spherical_hankel1`.

    This function describes a solution to the scalar Helmholtz wave
    equation in spherical coordinates.

    Args:
        l (int or array_like): Degree :math:`l \geq 0`.
        m (int or array_like): Order :math:`|m| \leq l`.
        kr (float or complex or array_like): Distance in units of the wavenumber in the medium :math:`kr`.
        theta (float or complex or array_like): Polar angle (rad).
        phi (float or array_like): Azimuthal angle (rad).

    Returns:
        complex
    """
    return sc.sph_harm(m, l, phi, theta) * sc.spherical_hankel1(l, kr)

def ssw_psi(l, m, x, y, z, theta, phi, k):
    r"""Far-field amplitude of singular scalar spherical wave :math:`\psi`

    It is defined as

    .. math::

        \psi_{lm}(x, y, z, \theta, \varphi)
        = \frac{\mathrm{(-i)}^{l+1}}{k} 
        Y_{lm}(\theta, \varphi)
        \mathrm{e}^{-\mathrm{i} k \left(x \sin \theta \cos \varphi + y \sin \theta \sin \varphi + z \cos \theta\right)}

    with :py:obj:`treams.special.sph_harm`.

    This function describes a part of the solution to the vector Helmholtz wave
    equation in spherical coordinates that does not depend on :math:`r` 
    in the limit of :math:`kr \to +\infty`.

    Args:
        l (int or array_like): Degree :math:`l \geq 0`.
        m (int or array_like): Order :math:`|m| \leq l`.
        x (float or complex or array_like): X-coordinate of the source.
        y (float or complex or array_like): Y-coordinate of the source. Has units of x.
        z (float or complex or array_like): Z-coordinate of the source. Has units of x.
        theta (float or complex or array_like): Polar angle (rad).
        phi (float or array_like): Azimuthal angle (rad).
        k (float or complex or array_like): Angular wavenumber in the medium. Has units of 1/x.

    Returns:
        complex
     """
    return (sc.sph_harm(m, l, phi, theta)
             * np.exp(-1j * ((k * x * np.cos(phi) + k * y * np.sin(phi)) * np.sin(theta) 
                             + k * z * np.cos(theta)))
             * np.power(-1j, l + 1)
             * 1 / k)

def ssw_rPsi(l, m, kr, theta, phi):
    r"""Regular scalar spherical wave :math:`\Psi`

    The regular scalar spherical wave is defined as

    .. math::

        \Psi_{lm}^{(1)}(kr, \theta, \varphi)
        = j_l(kr) Y_{lm}(\theta, \varphi)

    with :py:obj:`treams.special.sph_harm` and :py:obj:`scipy.special.spherical_jn`.

    This function describes a solution to the scalar Helmholtz wave
    equation in spherical coordinates.

    Args:
        l (int or array_like): Degree :math:`l \geq 0`.
        m (int or array_like): Order :math:`|m| \leq l`.
        kr (float or complex or array_like): Distance in units of the wavenumber in the medium :math:`kr`.
        theta (float or complex or array_like): Polar angle (rad).
        phi (float or array_like): Azimuthal angle (rad).

    Returns:
        complex
    """
    return sc.sph_harm(m, l, phi, theta) * sc.spherical_jn(l, kr)


def scw_Psi(kz, m, krr, phi, z):
    r"""Singular scalar cylindrical wave :math:`\Psi`

    The singular scalar cylindrical wave is defined as

    .. math::

        \Psi_{k_z m}^{(3)}(k_\rho \rho, \varphi, z)
        = H_m^{(1)}(k_\rho \rho) \mathrm e^{\mathrm i (k_z z + m \varphi)}

    using :py:obj:`scipy.special.hankel1`.

    This function describes a solution to the scalar Helmholtz wave
    equation in cylindrical coordinates.

    Args:
        kz (float or array_like): Z-component of the wave vector.
        m (int or array_like): Order :math:`m`.
        krr (float or complex or array_like): Radial distance in units of the radial component of the wave vector in the medium :math:`k_\rho \rho`.
        phi (float or array_like): Azimuthal angle (rad).
        z (float or array_like): Z-coordinate. Has units of 1/kz.

    Returns:
        complex
    """
    return np.exp(1j * (m * phi + kz * z)) * sc.hankel1(m, krr)

def scw_psi(kz, m, x, y, phi, z, krho):
    r"""Far-field amplitude of singular scalar cylindrical wave :math:`\psi`

    It is defined as

    .. math::

        \psi_{k_z m}(x, y, \varphi)
        = \sqrt{\frac{2}{\pi k_{\rho}}} \mathrm{i}^{-m} \mathrm{e}^{-\mathrm{i} \pi/4}
        \mathrm{e}^{-\mathrm{i} k_{\rho} \left(x \cos \varphi + y \sin \varphi\right)}
        \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}

    This function describes a part of the solution to the vector Helmholtz wave
    equation in cylindrical coordinates that does not depend on :math:`\rho` 
    in the limit of :math:`k_{\rho} \rho \to +\infty`.

    Args:
        kz (float or array_like): Z-component of the wave vector.
        m (int or array_like): Order.
        x (float or complex or array_like): X-coordinate of the source.
        y (float or complex or array_like): Y-coordinate of the source. Has units of x.
        phi (float or array_like): Azimuthal angle (rad).
        z (float or complex or array_like): Z-coordinate of the source. Has units of 1/kz.
        krho (float or complex or array_like): Radial component of the wave vector. Has units of 1/x.

    Returns:
        complex
     """
    m = np.asarray(m)
    krho = np.asarray(krho)
    with np.errstate(divide="ignore", invalid="ignore"):
        res = (
            np.exp(1j * m * phi + 1j * kz * z)
            * np.sqrt(2 / (np.pi * krho))
            * np.exp(-1j * krho * (x * np.cos(phi) + y * np.sin(phi)))
            * np.power(-1j, m)
            * np.exp(-1j * np.pi / 4)
        )
     # Evanescent modes (complex k_rho) do not contribute to the far field
    return np.where(krho.imag != 0, 0j, res)

def scw_rPsi(kz, m, krr, phi, z):
     r"""Regular scalar cylindrical wave :math:`\Psi`

     The regular scalar cylindrical wave is defined as

     .. math::

         \Psi_{k_z m}^{(1)}(k_\rho \rho, \varphi, z)
        = J_m(k_\rho \rho) \mathrm e^{\mathrm i (k_z z + m \varphi)}

     using :py:obj:`scipy.special.jv`.

     This function describes a solution to the scalar Helmholtz wave
     equation in cylindrical coordinates.

     Args:
         kz (float or array_like): Z-component of the wave vector.
         m (int or array_like): Order.
         krr (float or complex or array_like): Radial distance in units of the radial component of the wave vector in the medium :math:`k_\rho \rho`.
         phi (float or array_like): Azimuthal angle (rad).
         z (float or array_like): Z-coordinate. Has units of 1/kz.

     Returns:
         complex
     """
     return np.exp(1j * (m * phi + kz * z)) * sc.jv(m, krr)


def vsw_L(l, m, kr, theta, phi):
    r"""Singular longitudinal vector spherical wave L

    The singular vector spherical wave is defined as

    .. math::

        \mathbf{L}_{lm}^{(3)}(kr, \theta, \varphi)
        = -\mathrm{i} \left[{h_l^{(1)}}'(kr) \mathbf{Z}_{lm}(\theta, \varphi) 
        + \sqrt{l(l+1)} \frac{h_l^{(1)}(kr)}{kr} \mathbf{Y}_{lm}(\theta, \varphi) \right]

    with :py:obj:`treams.special.vsh_Z`, :py:obj:`treams.special.vsh_Y`, 
    :py:obj:`treams.special.spherical_hankel1`, and :py:obj:`treams.special.spherical_hankel1_d`.

    This function describes a longitudinal solution to the vector Helmholtz wave
    equation in spherical coordinates.

    Args:
        l (int or array_like): Degree :math:`l \geq 0`.
        m (int or array_like): Order :math:`|m| \leq l`.
        kr (float or complex or array_like): Distance in units of the wave number in the medium :math:`kr`.
        theta (float or complex or array_like): Polar angle (rad).
        phi (float or array_like): Azimuthal angle (rad).

    Returns:
        complex array: Components in spherical coordinates in
        the last dimension, which has length 3.
     """
    k_l = np.sqrt(l * (l + 1))[..., None]
    h = sc.spherical_hankel1(l, kr)
    res = (
        sc.vsh_Z(l, m, theta, phi) * sc.spherical_hankel1_d(l, kr)[..., None]
        + k_l * sc.vsh_Y(l, m, theta, phi) * (h / kr)[..., None]
    )
    return -1j * res

def vsw_rL(l, m, kr, theta, phi):
    r"""Regular longitudinal vector spherical wave L

    The regular vector spherical wave is defined as

    .. math::

        \mathbf{L}_{lm}^{(1)}(kr, \theta, \varphi)
        = -\mathrm{i} \left[j_l'(kr) \mathbf{Z}_{lm}(\theta, \varphi) 
        + \sqrt{l(l+1)}\frac{j_l(kr)}{kr} \mathbf{Y}_{lm}(\theta, \varphi) \right]

    with :py:obj:`treams.special.vsh_Z`, :py:obj:`treams.special.vsh_Y`, 
    and :py:obj:`scipy.special.spherical_jn`.

    This function describes a longitudinal solution to the vector Helmholtz wave
    equation in spherical coordinates.

    Args:
        l (int or array_like): Degree :math:`l \geq 0`.
        m (int or array_like): Order :math:`|m| \leq l`.
        kr (float or complex or array_like): Distance in units of the wave number in the medium :math:`kr`.
        theta (float or complex or array_like): Polar angle (rad).
        phi (float or array_like): Azimuthal angle (rad).

    Returns:
        complex array: Components in spherical coordinates in
        the last dimension, which has length 3.
     """
    l, m, kr, theta, phi = np.broadcast_arrays(l, m, kr, theta, phi)
    # At the origin the angles are undefined; the limit is evaluated at theta = phi = 0.
    at_origin = kr == 0
    theta = np.where(at_origin, 0, theta)
    phi = np.where(at_origin, 0, phi)
    # Limits of j_l'(x) and j_l(x) / x for x -> 0: 1/3 for l = 1, zero otherwise.
    lim = np.where(l == 1, 1 / 3, 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        jn_d = np.where(at_origin, lim, sc.spherical_jn_d(l, kr))
        jn_x = np.where(at_origin, lim, sc.spherical_jn(l, kr) / kr)
    res = (
        sc.vsh_Z(l, m, theta, phi) * jn_d[..., None]
        + np.sqrt(l * (l + 1))[..., None] * sc.vsh_Y(l, m, theta, phi) * jn_x[..., None]
    )
    return -1j * res


def vcw_L(kz, m, krr, phi, z, krho, k):
     r"""Singular longitudinal vector cylindrical wave L

     The singular vector cylindrical wave is defined as

     .. math::

         \mathbf{L}_{k_z m}^{(3)}(k_\rho \rho, \varphi, z)
        = \left[\frac{k_{\rho}}{k} {H^{(1)}_m(k_\rho \rho)}' \hat{\boldsymbol{\rho}} 
        + \mathrm{i}\frac{m k_{\rho}}{k}\frac{H^{(1)}_m(k_\rho \rho)}{k_{\rho}\rho} \hat{\boldsymbol{\varphi}} 
        + \frac{\mathrm{i} k_z}{k}H^{(1)}_m(k_\rho \rho) \hat{\mathbf{z}} \right] 
        \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}

     using :py:obj:`treams.special.hankel1_d` and :py:obj:`scipy.special.hankel1`.

     This function describes a longitudinal solution to the vector Helmholtz wave
     equation in cylindrical coordinates.

     Args:
         kz (float or array_like): Z-component of the wave vector.
         m (int or array_like): Order.
         krr (float or complex or array_like): Radial distance in units of the radial component of the wave vector in the medium :math:`k_\rho \rho`.
         phi (float or array_like): Azimuthal angle (rad).
         z (float or array_like): Z-coordinate. Has units of 1/kz.
         krho (float or complex or array_like): Radial component of the wave vector. Has units of kz.
         k (float or complex or array_like): Angular wavenumber in the medium. Has units of kz.

     Returns:
         complex array: Components in cylindrical coordinates in
         the last dimension, which has length 3.
      """
     h = sc.hankel1(m, krr)
     res = np.stack(
         np.broadcast_arrays(
             sc.hankel1_d(m, krr) * krho / k,
             1j * m * h / krr * krho / k,
             1j * kz / k * h,
         ),
         axis=-1,
     )
     return res * np.exp(1j * (m * phi + kz * z))[..., None]


def vcw_rL(kz, m, krr, phi, z, krho, k):
     r"""Regular longitudinal vector cylindrical wave L

     The regular vector cylindrical wave is defined by

     .. math::

         \mathbf{L}_{k_z m}^{(1)}(k_\rho \rho, \varphi, z)
        = \left[\frac{k_{\rho}}{k} {J_m(k_\rho \rho)}' \hat{\boldsymbol{\rho}} 
        + \mathrm{i}\frac{m k_{\rho}}{k}\frac{J_m(k_\rho \rho)}{k_{\rho}\rho} \hat{\boldsymbol{\varphi}} 
        + \frac{\mathrm{i} k_z}{k}J_m(k_\rho \rho) \hat{\mathbf{z}} \right] 
        \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}

     using :py:obj:`treams.special.jv_d` and :py:obj:`scipy.special.jv`.

     This function describes a longitudinal solution to the vector Helmholtz wave
     equation in cylindrical coordinates.

     Args:
         kz (float or array_like): Z-component of the wave vector.
         m (int or array_like): Order.
         krr (float or complex or array_like): Radial distance in units of the radial component of the wave vector in the medium :math:`k_\rho \rho`.
         phi (float or array_like): Azimuthal angle (rad).
         z (float or array_like): Z-coordinate. Has units of 1/kz.
         krho (float or complex or array_like): Radial component of the wave vector. Has units of kz.
         k (float or complex or array_like): Angular wavenumber in the medium. Has units of kz.

     Returns:
         complex array: Components in cylindrical coordinates in
         the last dimension, which has length 3.
     """
     m = np.asarray(m)
     on_axis = np.asarray(krr) == 0
     # Limits of J_m'(x) and J_m(x) / x for x -> 0: m / 2 for |m| = 1, zero otherwise.
     lim = np.where(np.abs(m) == 1, 0.5 * m, 0)
     with np.errstate(divide="ignore", invalid="ignore"):
         jv_d = np.where(on_axis, lim, sc.jv_d(m, krr))
         jv_x = np.where(on_axis, lim, sc.jv(m, krr) / krr)
     res = np.stack(
         np.broadcast_arrays(
             jv_d * krho / k,
             1j * m * jv_x * krho / k,
             1j * kz / k * sc.jv(m, krr),
         ),
         axis=-1,
     )
     return res * np.exp(1j * (m * phi + kz * z))[..., None]

def tl_ssw(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs):
    r"""tl_ssw(lambda_, mu, l, m, kr, theta, phi)
    
    Singular translation coefficient of scalar and longitudinal vector spherical waves.

    Defined by

    .. math::

        C_{\lambda \mu l m}^{(3)}(kr, \theta, \varphi) = (-1)^m \mathrm{i}^{\lambda-l} \sqrt{4\pi (2l+1)(2\lambda+1)} \\
        \cdot \sum_q \mathrm{i}^{q} \sqrt{2q+1} h^{(1)}_q(kr) Y_{q,m-\mu}(\theta,\varphi) \\
        \cdot \begin{pmatrix}
            l& \lambda & q \\
            m & -\mu & \mu-m
        \end{pmatrix}
        \begin{pmatrix}
            l& \lambda & q \\
            0& 0& 0
        \end{pmatrix}

    with the Wigner 3j-symbols (:py:obj:`treams.special.wigner3j`), the spherical Hankel
    functions (:py:obj:`treams.special.spherical_hankel1`), and the spherical harmonics
    (:py:obj:`treams.special.sph_harm`). The summation runs over all
    :math:`q \in \{\lambda + l, \lambda + l - 2, \dots, \max(|\lambda - l|, |\mu - m|)\}`.

    These coefficients are used to translate from scattered to incident modes.

    Args:
        lambda (int or array_like): Degree of the output mode.
        mu (int or array_like): Order of the output mode.
        l (int or array_like): Degree of the input mode.
        m (int or array_like): Order of the input mode.
        x (float or complex or array_like): Translation in units of the wavenumber in the medium.
        theta (float or complex or array_like): Polar angle (rad).
        phi (float or array_like): Azimuthal angle (rad).

    Returns:
        complex
    """
    pref = np.power(-1, np.abs(m)) * np.sqrt((2 * l + 1) * (2 * lambda_ + 1)) * cmath.exp(1j * (m - mu) * phi)
    res = 0.
    max_ = np.max([np.abs(int(lambda_) - int(l)), np.abs(int(m) - int(mu))])
    min_ = int(l) +  int(lambda_)
    for p in range(min_, max_ - 1, -2):
        res += (
            sc._tl_vsw_helper(l, m, lambda_, -mu, p, p)
            * sc.spherical_hankel1(p, kr)
            * sc.lpmv(m - mu, p, np.cos(theta), *args, **kwargs)
        )
    return res * pref


def tl_ssw_r(lambda_, mu, l, m, kr, theta, phi, *args, **kwargs):
    r"""tl_ssw_r(lambda_, mu, l, m, kr, theta, phi)
    
    Regular translation coefficient of scalar and longitudinal vector spherical waves.

    Defined by

    .. math::

        C_{\lambda \mu l m}^{(1)}(kr, \theta, \varphi) = (-1)^m \mathrm{i}^{\lambda-l} \sqrt{4\pi (2l+1)(2\lambda+1)} \\
        \cdot \sum_q \mathrm{i}^{q} \sqrt{2q+1} j_q(kr) Y_{q,m-\mu}(\theta,\varphi) \\
        \cdot \begin{pmatrix}
            l& \lambda & q \\
            m & -\mu & \mu-m
        \end{pmatrix}
        \begin{pmatrix}
            l& \lambda & q \\
            0& 0& 0
        \end{pmatrix}

    with the Wigner 3j-symbols (:py:obj:`treams.special.wigner3j`), the spherical Bessel
    functions (:py:obj:`scipy.special.spherical_jn`), and the spherical harmonics
    (:py:obj:`treams.special.sph_harm`). The summation runs over all
    :math:`q \in \{\lambda + l, \lambda + l - 2, \dots, \max(|\lambda - l|, |\mu - m|)\}`.

    These coefficients are used to translate from scattered to scattered modes and
    from incident to incident modes.

    Args:
        lambda (int or array_like): Degree of the output mode.
        mu (int or array_like): Order of the output mode.
        l (int or array_like): Degree of the input mode.
        m (int or array_like): Order of the input mode.
        kr (float or complex or array_like): Translation in units of the wavenumber in the medium.
        theta (float or complex or array_like): Polar angle (rad).
        phi (float or array_like): Azimuthal angle (rad).

    Returns:
        complex
    """
    pref = np.power(-1, np.abs(m)) * np.sqrt((2 * l + 1) * (2 * lambda_ + 1)) * cmath.exp(1j * (m - mu) * phi)
    res = 0.
    max_ = np.max([np.abs(int(lambda_) - int(l)), np.abs(int(m) - int(mu))])
    min_ = int(l) + int(lambda_)
    for p in range(min_, max_ - 1, -2):
        res += (
            sc._tl_vsw_helper(l, m, lambda_, -mu, p, p)
            * sc.spherical_jn(p, kr)
            * sc.lpmv(m - mu, p, np.cos(theta), *args, **kwargs)
        )
    return res * pref

def tl_scw(kz1, mu, kz2, m, krr, phi, z, *args, **kwargs):
    r"""tl_scw(kz1, mu, kz2, m, krr, phi, z)
    
    Singular translation coefficient of scalar and longitudinal vector cylindrical waves.

    Defined by

    .. math::

        \begin{cases}
            H_{m - \mu}^{(1)}(k_\rho \rho) \mathrm e^{\mathrm i ((m - \mu) \varphi + k_{z,1}z)}
            & \text{if }k_{z,1} = k_{z,2} \\
            0 & \text{otherwise}
        \end{cases}

    where :math:`H_{m - \mu}^{(1)}` are the Hankel functions (:py:obj:`scipy.special.hankel1`).

    These coefficients are used to translate from scattered to incident modes.

    Args:
        kz1 (float or array_like): Z-component of the wave vector of the output mode.
        mu (int or array_like): Order of the output mode
        kz2 (float or array_like): Z-component of the wave vector of the input mode. Has units of kz1.
        m (int or array_like): Order of the input mode.
        krr (float or complex or array_like): Translation in radial direction in units of 
            the radial component of the wave vector in the medium :math:`k_\rho \rho`.
        phi (float or array_like): Azimuthal angle.
        z (float or array_like): Translation in the z-direction. Has units of 1/kz1.

    Returns:
        complex
    """
    if kz1 != kz2:
        return 0.+0.j
    return sc.hankel1(m - mu, krr) * np.exp(1j * ((m - mu) * phi + kz1 * z), *args, **kwargs)


def tl_scw_r(kz1, mu, kz2, m, krr, phi, z, *args, **kwargs):
    r"""tl_scw_r(kz1, mu, kz2, m, krr, phi, z)
    
    Regular translation coefficient of scalar and longitudinal vector cylindrical waves.

    Defined by

    .. math::

        \begin{cases}
            J_{m - \mu}(k_\rho \rho) \mathrm e^{\mathrm i ((m - \mu) \varphi + k_{z,1}z)}
            & \text{if }k_{z,1} = k_{z,2} \\
            0 & \text{otherwise}
        \end{cases}

    where :math:`J_{m - \mu}` are the Bessel functions (:py:obj:`scipy.special.jv`).

    These coefficients are used to translate from scattered to scattered modes or
    from incident to incident modes.

    Args:
        kz1 (float or array_like): Z-component of the wave vector of the output mode.
        mu (int or array_like): Order of the output mode
        kz2 (float or array_like): Z-component of the wave vector of the input mode. Has units of kz1.
        m (int or array_like): Order of the input mode.
        krr (float or complex or array_like): Translation in radial direction in units of 
            the radial component of the wave vector in the medium :math:`k_\rho \rho`.
        phi (float or array_like): Azimuthal angle.
        z (float or array_like): Translation in the z-direction. Has units of 1/kz1.

    Returns:
        complex
    """
    if kz1 != kz2:
        return 0.+0.j
    return sc.jv(m - mu, krr) * np.exp(1j * ((m - mu) * phi + kz1 * z), *args, **kwargs)
