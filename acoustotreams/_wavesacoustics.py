import treams.special as sc
import numpy as np
import cmath
import scipy.special as ss

def spw_Psi(kx, ky, kz, x, y, z):
    r"""Scalar plane wave :math:`\Psi` 

    The scalar plane wave is defined by

    .. math::

        \psi_{\mathbf k}(\mathbf r)
        = \mathrm e^{\mathrm i \mathbf k \mathbf r}.

    This function describes a solution to the scalar Helmholtz wave
    equation.

    Args:
        kx (float or complex, array_like): X component of the wave vector
        ky (float or complex, array_like): Y component of the wave vector
        kz (float or complex, array_like): Z component of the wave vector
        x (float, array_like): Y coordinate
        y (float, array_like): X coordinate
        z (float, array_like): Z coordinate

    Returns:
        complex
    
    """
    return np.exp(1j * (kx * x + ky * y + kz * z))


def vpw_L(kx, ky, kz, x, y, z):
    r"""Longitudinal vector plane wave L
    
    The longitudinal vector plane wave is defined by

    .. math::

       \mathbf L_{\mathbf k}(\mathbf r)
        = \frac{\mathbf k}{k} \mathrm e^{\mathrm i \mathbf k \mathbf r}.

    This function describes a longitudinal solution to the vector Helmholtz wave
    equation.

    Args:
        kx (float or complex, array_like): X component of the wave vector
        ky (float or complex, array_like): Y component of the wave vector
        kz (float or complex, array_like): Z component of the wave vector
        x (float, array_like): Y coordinate
        y (float, array_like): X coordinate
        z (float, array_like): Z coordinate

    Returns:
        complex, 3-array
    
    """
    k = np.sqrt(kx * kx + ky * ky + kz * kz)
    phase = np.exp(1j * (kx * x + ky * y + kz * z))
    if k == 0:
        return np.asarray([np.nan, np.nan, np.nan]).T
    else:
        return 1j * np.asarray([kx * phase, ky * phase, kz * phase], np.complex128).T / k

def ssw_Psi(l, m, kr, theta, phi):
    r"""Singular scalar spherical wave :math:`\Psi`

    The singular scalar spherical wave is defined by

    .. math::

        \Psi_{lm}^{(3)}(x, \theta, \varphi)
        = h_l^{(1)}(x) Y_{lm}(\theta, \varphi)

    with :func:`acoustotreams.sph_harm`, and :func:`acoustotreams.spherical_hankel1`.

    This function describes a solution to the scalar Helmholtz wave
    equation.

    Args:
        l (int, array_like): Degree :math:`l \geq 0`
        m (int, array_like): Order :math:`|m| \leq l`
        x (float or complex, array_like): Distance in units of the wave number in the medium :math:`kr`
        theta (float or complex, array_like): Polar angle
        phi (float, array_like): Azimuthal angle

    Returns:
        complex
    """
    return sc.sph_harm(m, l, phi, theta) * sc.spherical_hankel1(l, kr)

def ssw_psi(l, m, x, y, z, theta, phi, k):
     r"""Far-field amplitude of singular scalar spherical wave :math:`\psi`
     Defined by

    .. math::

        \psi_{lm}(x, y, z, \theta, \varphi)
        = \frac{\mathrm{i}^{-l-1}}{k} 
        Y_{lm}(\theta, \varphi)
        \mathrm{e}^{-\mathrm{i} k \left(x \sin \theta \cos \varphi + y \sin \theta \sin \varphi + z \cos \theta\right)}

    with :func:`acoustotreams.sph_harm`.

    This function describes a solution to the vector Helmholtz wave
    equation in the limit of :math:`kr \to +\infty`.

    Args:
        l (int, array_like): Degree :math:`l \geq 0`
        m (int, array_like): Order :math:`|m| \leq l`
        x (float or complex, array_like): X coordinate of the source
        y (float or complex, array_like): Y coordinate of the source
        z (float or complex, array_like): Z coordinate of the source
        theta (float or complex, array_like): Polar angle
        phi (float, array_like): Azimuthal angle
        k (float or complex, array_like): Wave number in the medium

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

     The regular scalar spherical wave is defined by

    .. math::

        \Psi_{lm}^{(1)}(x, \theta, \varphi)
        = j_l(x) Y_{lm}(\theta, \varphi)

    with :func:`acoustotreams.sph_harm`, and :func:`acoustotreams.spherical_jn`.

    This function describes a solution to the scalar Helmholtz wave
    equation.

    Args:
        l (int, array_like): Degree :math:`l \geq 0`
        m (int, array_like): Order :math:`|m| \leq l`
        x (float or complex, array_like): Distance in units of the wave number in the medium :math:`kr`
        theta (float or complex, array_like): Polar angle
        phi (float, array_like): Azimuthal angle

    Returns:
        complex
    """
     return sc.sph_harm(m, l, phi, theta) * sc.spherical_jn(l, kr)


def scw_Psi(kz, m, krr, phi, z):
     r"""Singular scalar cylindrical wave :math:`\Psi`

     The singular scalar cylindrical wave is defined by

     .. math::

         \psi_{k_z m}^{(3)}(x_\rho, \varphi, z)
        = H_m^{(1)}(x_\rho) \mathrm e^{\mathrm i (k_z z + m \varphi)}

     using :func:`acoustotreams.hankel1`.

     This function describes a solution to the scalar Helmholtz wave
     equation.

     Args:
         kz (float, array_like): Z component of the wave vector
         m (int, array_like): Order :math:`|m| \leq l`
         xrho (float or complex, array_like): Radial distance in units of the radial component of the wave number in the medium :math:`k_\rho \rho`
         phi (float, array_like): Azimuthal angle
         z (float, array_like): Z coordinate

     Returns:
         complex
      """
     return np.exp(1j * (m * phi + kz * z)) * sc.hankel1(m, krr)

def scw_psi(kz, m, x, y, phi, z, krho):
     r"""Far-field amplitude of singular scalar cylindrical wave :math:`\psi`

     Defined by

    .. math::

        \psi_{k_z m}(x, y, \varphi)
        = \sqrt{\frac{2}{\pi k_{\rho}}} \mathrm{i}^{-m} \mathrm{e}^{-\mathrm{i} \pi/4}
        \mathrm{e}^{-\mathrm{i} k_{\rho} \left(x \cos \varphi + y \sin \varphi\right)}
        \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}

    This function describes a solution to the vector Helmholtz wave
    equation in the limit of :math:`k_{\rho} \rho \to +\infty`.

    Args:
        kz (float, array_like): Z component of the wave vector
        m (int, array_like): Order
        x (float or complex, array_like): X coordinate of the source
        y (float or complex, array_like): Y coordinate of the source
        phi (float, array_like): Azimuthal angle
        z (float or complex, array_like): Z coordinate of the source
        krho (float or complex, array_like): Radial component of the wave vector

    Returns:
        complex
     """
     return (np.exp(1j * m * phi + 1j * kz * z) * 
             np.sqrt(2 / (np.pi * krho)) *
             np.exp(-1j * ((krho * x * np.cos(phi) + krho * y * np.sin(phi)))) *
             np.power(-1j, m) * 
             np.exp(-1j * np.pi/4))

def scw_rPsi(kz, m, krr, phi, z):
     r"""Regular scalar cylindrical wave :math:`\Psi`

     The regular scalar cylindrical wave is defined by

     .. math::

         \psi_{k_z m}^{(1)}(x_\rho, \varphi, z)
        = J_m(x_\rho) \mathrm e^{\mathrm i (k_z z + m \varphi)}

     using :func:`acoustotreams.jv`.

     This function describes a solution to the scalar Helmholtz wave
     equation.

     Args:
         kz (float, array_like): Z component of the wave vector
         m (int, array_like): Order :math:`|m| \leq l`
         xrho (float or complex, array_like): Radial distance in units of the radial component of the wave number in the medium :math:`k_\rho \rho`
         phi (float, array_like): Azimuthal angle
         z (float, array_like): Z coordinate

     Returns:
         complex
     """
     return np.exp(1j * (m * phi + kz * z)) * sc.jv(m, krr)


def vsw_L(l, m, kr, theta, phi):
     r"""Singular longitudinal vector spherical wave L

     The singular vector spherical wave is defined by

    .. math::

        \mathbf{L}_{lm}^{(3)}(x, \theta, \varphi)
        = -\mathrm{i} \left[{h_l^{(1)}}'(x) \mathbf{Z}_{lm}(\theta, \varphi) 
        + \sqrt{l(l+1)} \frac{h_l^{(1)}(x)}{x} \mathbf{Y}_{lm}(\theta, \varphi) \right]

    with :func:`acoustotreams.vsh_Z`, :func:`acoustotreams.vsh_Y`, 
    :func:`acoustotreams.spherical_hankel1`, and :func:`acoustotreams.spherical_hankel1`.

    This function describes a solution to the vector Helmholtz wave
    equation.

    Args:
        l (int, array_like): Degree :math:`l \geq 0`
        m (int, array_like): Order :math:`|m| \leq l`
        x (float or complex, array_like): Distance in units of the wave number in the medium :math:`kr`
        theta (float or complex, array_like): Polar angle
        phi (float, array_like): Azimuthal angle

    Returns:
        complex, 3-array
     """
     return np.transpose(sc.vsh_Z(l, m, theta, phi).T * sc.spherical_hankel1_d(l, kr)  + np.sqrt(
         l * (l + 1)
     ) * sc.vsh_Y(l, m, theta, phi).T * sc.spherical_hankel1(l, kr) / kr) * (-1.0j)


def vsw_l(l, m, x, y, z, theta, phi, k):
     r"""Far-field amplitude of singular longitudinal vector spherical wave l

     Defined by

    .. math::

        \mathbf{l}_{lm}(x, y, z, \theta, \varphi)
        = \frac{\mathrm{i}^{-l}}{k} 
        Y_{lm}(\theta, \varphi)\hat{\mathbf{r}}
        \mathrm{e}^{-\mathrm{i} k \left(x \sin \theta \cos \varphi + y \sin \theta \sin \varphi + z \cos \theta\right)}

    with :func:`acoustotreams.sph_harm`.

    This function describes a solution to the vector Helmholtz wave
    equation in the limit of :math:`kr \to +\infty`.

    Args:
        l (int, array_like): Degree :math:`l \geq 0`
        m (int, array_like): Order :math:`|m| \leq l`
        x (float or complex, array_like): X coordinate of the source
        y (float or complex, array_like): Y coordinate of the source
        z (float or complex, array_like): Z coordinate of the source
        theta (float or complex, array_like): Polar angle
        phi (float, array_like): Azimuthal angle
        k (float or complex, array_like): Wave number in the medium

    Returns:
        complex, 3-array
     """
     return np.transpose(
         np.array(
             [
                 sc.sph_harm(m, l, phi, theta)
                 * np.exp(-1j * ((k * x * np.cos(phi) + k * y * np.sin(phi)) * np.sin(theta) 
                                 + k * z * np.cos(theta)))
                 * np.power(-1j, l)
                 * 1 / k, 
                 np.zeros(len(l)).T,
                 np.zeros(len(l)).T
            ]
        )
    )


def vsw_rL(l, m, kr, theta, phi):
     r"""Regular longitudinal vector spherical wave L

     The regular vector spherical wave is defined by

    .. math::

        \mathbf{L}_{lm}^{(1)}(x, \theta, \varphi)
        = -\mathrm{i} \left[j_l'(x) \mathbf{Z}_{lm}(\theta, \varphi) 
        + \sqrt{l(l+1)}\frac{j_l(x)}{x} \mathbf{Y}_{lm}(\theta, \varphi) \right]

    with :func:`acoustotreams.vsh_Z`, :func:`acoustotreams.vsh_Y`, 
    and :func:`acoustotreams.spherical_jn`.

    This function describes a solution to the vector Helmholtz wave
    equation.

    Args:
        l (int, array_like): Degree :math:`l \geq 0`
        m (int, array_like): Order :math:`|m| \leq l`
        x (float or complex, array_like): Distance in units of the wave number in the medium :math:`kr`
        theta (float or complex, array_like): Polar angle
        phi (float, array_like): Azimuthal angle

    Returns:
        complex, 3-array
     """
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
     r"""Singular longitudinal vector cylindrical wave L

     The singular vector cylindrical wave is defined by

     .. math::

         \mathbf{L}_{k_z m}^{(3)}(x_\rho, \varphi, z)
        = \left[\frac{k_{\rho}}{k} {H^{(1)}_m(x_\rho)}' \hat{\boldsymbol{\rho}} 
        + \mathrm{i}\frac{m k_{\rho}}{k}\frac{H^{(1)}_m(x_\rho)}{k_{\rho}\rho} \hat{\boldsymbol{\varphi}} 
        + \frac{\mathrm{i} k_z}{k}H^{(1)}_m(x_\rho) \hat{\mathbf{z}} \right] 
        \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}

     using :func:`acoustotreams.hankel1_d` and :func:`acoustotreams.hankel1`.

     This function describes a longitudinal solution to the vector Helmholtz wave
     equation.

     Args:
         kz (float, array_like): Z component of the wave vector
         m (int, array_like): Order
         xrho (float or complex, array_like): Radial distance in units of the radial component of the wave number in the medium :math:`k_\rho \rho`
         phi (float, array_like): Azimuthal angle
         z (float, array_like): Z coordinate
         krho (float or complex, array_like): Radial component of the wave vector
         k (float or complex, array_like): Wavenumber in the medium 

     Returns:
         complex, 3-array
      """
     return np.transpose(
         np.array(
             [
                 sc.hankel1_d(m, krr) * krho / k  * np.exp(1j * (m * phi + kz * z)),
                 1j * m * sc.hankel1(m, krr) / krr * krho / k  * np.exp(1j * (m * phi + kz * z)),
                 1j * kz / k * sc.hankel1(m, krr) * np.exp(1j * (m * phi + kz * z)),
             ]
         )
     )

def vcw_l(kz, m, x, y, phi, z, krho, k):
     r"""Far-field amplitude of longitudinal singular vector cylindrical wave l

     Defined by

    .. math::

        \mathbf{l}_{k_z m}(x, y, \varphi, z,)
        = \sqrt{\frac{2}{\pi k_{\rho}}} \mathrm{i}^{1-m} \mathrm{e}^{-\mathrm{i} \pi/4}
        \left(\frac{k_{\rho}}{k} \hat{\boldsymbol{\rho}} 
         + \frac{k_{z}}{k} \hat{\mathbf{z}} \right)
        \mathrm{e}^{-\mathrm{i} k_{\rho} \left(x \cos \varphi + y \sin \varphi\right)}
        \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}

    This function describes a solution to the vector Helmholtz wave
    equation in the limit of :math:`k_{\rho} \rho \to +\infty`.

    Args:
        kz (float, array_like): Z component of the wave vector
        m (int, array_like): Order
        x (float or complex, array_like): X coordinate of the source
        y (float or complex, array_like): Y coordinate of the source
        phi (float, array_like): Azimuthal angle
        z (float or complex, array_like): Z coordinate of the source
        krho (float or complex, array_like): Radial component of the wave vector
        k (float or complex, array_like): Wave number in the medium

    Returns:
        complex, 3-array
     """
     return np.transpose(
         np.array(
             [
                np.sqrt(2 / (np.pi * krho)) *
                np.exp(-1j * ((krho * x * np.cos(phi) + krho * y * np.sin(phi)))) *
                np.power(-1j, m) * 
                np.exp(1j * (m * phi + kz * z)) *
                np.exp(-1j * np.pi/4) *  
                1j * krho / k,
                np.zeros(len(krho)).T,
                np.sqrt(2 / (np.pi * krho)) *
                np.exp(-1j * ((krho * x * np.cos(phi) + krho * y * np.sin(phi)))) *
                np.power(-1j, m) * 
                np.exp(1j * (m * phi + kz * z)) *
                np.exp(-1j * np.pi/4) *  
                1j * kz / k,
             ]
         )
     )

def vcw_rL(kz, m, krr, phi, z, krho, k):
     r"""Regular longitudinal vector cylindrical wave L

     The regular vector cylindrical wave is defined by

     .. math::

         \mathbf{L}_{k_z m}^{(1)}(x_\rho, \varphi, z)
        = \left[\frac{k_{\rho}}{k} {J_m(x_\rho)}' \hat{\boldsymbol{\rho}} 
        + \mathrm{i}\frac{m k_{\rho}}{k}\frac{J_m(x_\rho)}{k_{\rho}\rho} \hat{\boldsymbol{\varphi}} 
        + \frac{\mathrm{i} k_z}{k}J_m(x_\rho) \hat{\mathbf{z}} \right] 
        \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}

     using :func:`acoustotreams.jv_d` and :func:`acoustotreams.jv`.

     This function describes a longitudinal solution to the vector Helmholtz wave
     equation.

     Args:
         kz (float, array_like): Z component of the wave vector
         m (int, array_like): Order
         xrho (float or complex, array_like): Radial distance in units of the radial component of the wave number in the medium :math:`k_\rho \rho`
         phi (float, array_like): Azimuthal angle
         z (float, array_like): Z coordinate
         krho (float or complex, array_like): Radial component of the wave vector
         k (float or complex): Wavenumber in the medium 

     Returns:
         complex, 3-array
     """
     # todo: implement for k as array_like
     res = []
     for i, x in enumerate(krr):
            if x == 0:
                if (abs(m[i]) != 1) and (m[i] != 0):
                    val = np.zeros(3, complex).T
                    res.append(val)
                elif abs(m[i]) == 1:
                   val = np.array(
                                  [
                       0.5 * krho[i] / k,
                       0.5 * 1j * m[i] * krho[i] / k,
                       0. 
                       ]).T * np.exp(1j * (m[i] * phi[i] + kz[i] * z[i]))
                   res.append(val)
                elif m[i] == 0:
                    val = np.array(
                        [
                            0., 
                            0., 
                            1j * kz[i] / k * np.exp(1j * (m[i] * phi[i] + kz[i] * z[i]))
                            ]
                            ).T
                    res.append(val)
            else:      
                val = np.array(
                [
                    sc.jv_d(m[i], krr[i]) * krho[i] / k,
                    1j * m[i] * sc.jv(m[i], krr[i]) / krr[i] * krho[i] / k,
                    1j * kz[i] / k * sc.jv(m[i], krr[i]),
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
    r"""tl_ssw(lambda_, mu, l, m, kr, theta, phi)
    
    Singular translation coefficient of scalar and longitudinal vector spherical waves

    Definded by

    .. math::

        \alpha_{\lambda \mu l m}^{(3)}(x, \theta, \varphi) = (-1)^m \mathrm{i}^{\lambda-l} \sqrt{4\pi (2l+1)(2\lambda+1)} \\
        \cdot \sum_q \mathrm{i}^{q} \sqrt{2q+1} h^{(1)}_q(x) Y_{q,m-\mu}(\theta,\varphi) \\
        \cdot \begin{pmatrix}
            l& \lambda & q \\
            m & -\mu & \mu-m
        \end{pmatrix}
        \begin{pmatrix}
            l& \lambda & q \\
            0& 0& 0
        \end{pmatrix}

    with the Wigner 3j-symbols (:func:`acoustotreams.wigner3j`), the spherical Hankel
    functions (:func:`acoustotreams.spherical_hankel1`), and the spherical harmonics
    (:func:`acoustotreams.sph_harm`). The summation runs over all
    :math:`q \in \{\lambda + l, \lambda + l - 2, \dots, \max(|\lambda - l|, |\mu - m|)\}`.

    These coefficients are used to translate from scattered to incident modes.

    Args:
        lambda (integer, array_like): Degree of the destination mode
        mu (integer, array_like): Order of the destination mode
        l (integer, array_like): Degree of the source mode
        m (integer, array_like): Order of the source mode
        x (complex, array_like): Translation in units of the wave number in the medium
        theta (float or complex, array_like): Polar angle
        phi (float, array_like): Azimuthal angel

    Returns:
        complex
    """
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
    r"""tl_ssw_r(lambda_, mu, l, m, kr, theta, phi)
    
    Regular translation coefficient of scalar and longitudinal vector spherical waves

    Definded by

    .. math::

        \alpha_{\lambda \mu l m}^{(1)}(x, \theta, \varphi) = (-1)^m \mathrm{i}^{\lambda-l} \sqrt{4\pi (2l+1)(2\lambda+1)} \\
        \cdot \sum_q \mathrm{i}^{q} \sqrt{2q+1} j_q(x) Y_{q,m-\mu}(\theta,\varphi) \\
        \cdot \begin{pmatrix}
            l& \lambda & q \\
            m & -\mu & \mu-m
        \end{pmatrix}
        \begin{pmatrix}
            l& \lambda & q \\
            0& 0& 0
        \end{pmatrix}

    with the Wigner 3j-symbols (:func:`acoustotreams.wigner3j`), the spherical Bessel
    functions (:func:`acoustotreams.spherical_jn`), and the spherical harmonics
    (:func:`acoustotreams.sph_harm`). The summation runs over all
    :math:`q \in \{\lambda + l, \lambda + l - 2, \dots, \max(|\lambda - l|, |\mu - m|)\}`.

    These coefficients are used to translate from scattered to scattered modes and
    from incident to incident modes.

    Args:
        lambda (integer, array_like): Degree of the destination mode
        mu (integer, array_like): Order of the destination mode
        l (integer, array_like): Degree of the source mode
        m (integer, array_like): Order of the source mode
        x (complex, array_like): Translation in units of the wave number in the medium
        theta (float or complex, array_like): Polar angle
        phi (float, array_like): Azimuthal angel

    Returns:
        complex
    """
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
    r"""tl_scw(kz1, mu, kz2, m, krr, phi, z)
    
    Singular translation coefficient of scalar and longitudinal vector cylindrical waves

    Definded by

    .. math::

        \begin{cases}
            H_{m - \mu}^{(1)}(x_\rho) \mathrm e^{\mathrm i ((m - \mu) \varphi + k_{z,1}z)}
            & \text{if }k_{z,1} = k_{z,2} \\
            0 & \text{otherwise}
        \end{cases}

    where :math:`H_{m - \mu}^{(1)}` are the Hankel functions (:func:`acoustotreams.hankel1`).

    These coefficients are used to translate from scattered to scattered modes and
    from incident to incident modes.

    Args:
        kz1 (float, array_like): Z component of the destination mode's wave vector
        mu (integer, array_like): Order of the destination mode
        kz2 (float, array_like): Z component of the source mode's wave vector
        m (integer, array_like): Order of the source mode
        xrho (complex, array_like): Translation in radial direction in units of the radial component of the wave number in the medium :math:`k_\rho \rho`
        phi (float, array_like): Azimuthal angel
        z (float, array_like): Translation in z direction

    Returns:
        complex
    """
    if kz1 != kz2:
        return 0.+0.j
    return sc.hankel1(m - mu, krr) * np.exp(1j * ((m - mu) * phi + kz1 * z), *args, **kwargs)


def tl_scw_r(kz1, mu, kz2, m, krr, phi, z, *args, **kwargs):
    r"""tl_scw_r(kz1, mu, kz2, m, krr, phi, z)
    
    Regular translation coefficient of scalar and longitudinal vector cylindrical waves

    Definded by

    .. math::

        \begin{cases}
            J_{m - \mu}(x_\rho) \mathrm e^{\mathrm i ((m - \mu) \varphi + k_{z,1}z)}
            & \text{if }k_{z,1} = k_{z,2} \\
            0 & \text{otherwise}
        \end{cases}

    where :math:`J_{m - \mu}` are the Bessel functions (:func:`acoustotreams.jv`).

    These coefficients are used to translate from scattered to incident modes.

    Args:
        kz1 (float, array_like): Z component of the destination mode's wave vector
        mu (integer, array_like): Order of the destination mode
        kz2 (float, array_like): Z component of the source mode's wave vector
        m (integer, array_like): Order of the source mode
        xrho (complex, array_like): Translation in radial direction in units of the radial component of the wave number in the medium :math:`k_\rho \rho`
        phi (float, array_like): Azimuthal angel
        z (float, array_like): Translation in z direction

    Returns:
        complex
    """
    if kz1 != kz2:
        return 0.+0.j
    return sc.jv(m - mu, krr) * np.exp(1j * ((m - mu) * phi + kz1 * z), *args, **kwargs)