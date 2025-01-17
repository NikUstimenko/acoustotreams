======
Theory
======

These sections give an overview over theory that underlies the different aspects of acoustotreams.
For more in-depth information, there is also a list of associated :ref:`about:Publications`.

Master equations of linear monochromatic acoustics
==================================================

The master equations that describe propagation of pressure waves in nonelastic media
can be written  in frequency domain in the absence of external sources and shear (transverse) 
waves as [1]_

.. math::

    \begin{pmatrix}
    \boldsymbol{\nabla}& 0 \\
    0& \boldsymbol{\nabla} \cdot
    \end{pmatrix}
    \begin{pmatrix}
    p\\
    \mathbf{v}
    \end{pmatrix} = 
    \mathrm i \omega
    \begin{pmatrix}
    0&  \rho\\
    \frac{1}{\rho c^2}& 0
    \end{pmatrix}
    \begin{pmatrix}
    p\\
    \mathbf{v}
    \end{pmatrix}

where :math:`p = p(\mathbf{r}, \omega)` and :math:`\mathbf v = \mathbf v(mathbf r, \omega)`
are the pressure and velocity fields (:func:`acoustotreams.pfield`, :func:`acoustotreams.vfield`).
All these quantities are complex valued fields, that depend on the angular frequency :math:`\omega` and 
the position :math:`\mathbf r`. The material parameters are the speed of sound :math:`c`, and
the mass density :math:`\rho`. Conventionally, within acoustotreams the (air) wave number
:math:`k_0 = \frac{\omega}{c_0}` is generally used to express the frequency. Here :math:`c_0` is
the speed of sound in air being 343 m/s.

For the transformation to the time domain we use for a general function
:math:`f(\omega)`

.. math::

    f(t) = \int_{-\infty}^\infty \mathrm d t f(\omega) \mathrm e^{-\mathrm i \omega t}

as Fourier transformation convention, and thus the inverse transformation is

.. math::

    f(\omega)
    = \int_{-\infty}^\infty \frac{\mathrm d \omega}{2 \pi}
    f(t) \mathrm e^{\mathrm i \omega t}

Equation of motion for elastic objects
--------------------------------------

In acoustotreams we can also find the pressure and velocity fields in the presence of elastic objects
that support both pressure and shear waves while the background media support *no* shear waves. 
In this case we explicitly solve the equation for the displacement field :math:`\mathbf u = \frac{\mathrm i}{\omega}\mathbf v` [2]_

.. math::

    \omega^2 \mathbf u 
    + c_t^2 \boldsymbol \Delta \mathbf u 
    + (c^2 - c_t^2) \boldsymbol \nabla (\boldsymbol \nabla \cdot \mathbf u) 
    = 0

where :math:`\rho` is the mass density, :math:`c` is the speed of longitudinal pressure waves, 
and :math:`c_t` is the speed of transverse shear waves (:func:`acoustotreams.AcousticMaterial`).
Due to the requirement of isotropy these quantities are all scalar. In a medium with :math:`c_t \equiv 0`
the last equation coincides with the Helmholtz equations below and we can use the scalar description
of the waves.

Solutions to the scalar Helmholtz equation
==========================================

Instead of immediatly solving equations from above, we will study the
Helmholtz equation which is commonly encountered when studying wave phenomena first.
This section mainly relies on [3]_.

The scalar Helmholtz equation is

.. math::

    \left(\Delta + k^2 \right) p 
    = \boldsymbol{\nabla} \cdot (\boldsymbol{\nabla} p) - k^2 p
    = 0

where :math:`\Delta` is the Laplace operator. Note, that by applying the divergence operator
on equation :math:`\boldsymbol{\nabla}p=\mathrm i \omega \mathbf{v}` and using the equation
:math:`\boldsymbol{\nabla} \cdot \mathbf{v} = \frac{\mathrm i \omega}{\rho c^2}p`, the scalar
Helmholtz equation can be obtained. Alternatively, we can apply the gradient to the second equation
and get the vector Helmholtz equation

.. math::

    \left(\boldsymbol \Delta + k^2 \right) \mathbf{v} = \boldsymbol{\nabla} (\boldsymbol{\nabla} \cdot \mathbf v)
    - \boldsymbol{\nabla} \times \boldsymbol{\nabla} \times \mathbf{v}
    + k^2 \mathbf v
    = 0

If we denote :math:`\psi` as a solution to the scalar Helmholtz equation, solutions to the vector Helmholtz equation
can be constructed as follows

.. math::

    \mathbf L = k^{-1} \boldsymbol{\nabla} \psi \\
    \mathbf M = \boldsymbol{\nabla} \times (\mathbf c \psi) \\
    \mathbf N = \boldsymbol{\nabla} \times k^{-1} \boldsymbol{\nabla} \times (\mathbf c \psi)

where :math:`\mathbf c` is a steering vector that depends on the coordinate system
used for the solution :math:`\psi`. We will focus the following discussion on the three
cases of planar, cylindrical, and spherical solutions, where the coordinate systems are
chosen to be Cartesian, cylindrical, and spherical. 

The first type of solution is longitudinal waves that obey constraint :math:`\boldsymbol{\nabla} \times \mathbf L = 0`,
while the second and third ones are transverse waves with the constraint :math:`\boldsymbol{\nabla} \cdot \{\mathbf M,\mathbf N\} = 0`.
In the following, we will limit the discussion of the transverse waves, because acoustic waves in a medium with :math:`c_t \equiv 0` can
be described by scalar pressure fields and corresponding longitudinal fields.

Plane waves
-----------

In Cartesian coordinates the solution to the scalar Helmholtz equation are simple
plane waves 

.. math::

    \psi_{\mathbf k}(k, \mathbf r) = \mathrm e^{\mathrm i \mathbf k \mathbf r} 

given by :func:`acoustotreams.ssw_Psi` where the wave vector fulfils 
:math:`|\mathbf k|^2 = k_x^2 + k_y^2 + k_z^2 = k^2`. The corresponding 
longitudinal vector plane wave is

.. math::

    \mathbf L_{\mathbf k}(k, \mathbf r) 
    = \frac{k_x \mathbf{\hat x} + k_y \mathbf{\hat y} + k_z \mathbf{\hat z}}{k}
    = \mathbf{\hat r}_{\mathbf k}
    \mathrm e^{\mathrm i \mathbf k \mathbf r}

given by :func:`acoustotreams.vpw_L`. We normalized this wave by :math:`k` in the medium
such that it has unit strength for real-valued wave vectors.

Cylindrical waves
-----------------

The cylindrical solutions can be constructed mostly analogously to the plane waves.
The solutions in cylindrical coordinates are 

.. math::
    
    \psi^{(n)}_{k_z, m}(k, \mathbf r) 
    = 
    Z_m^{(n)}(k_\rho \rho) \mathrm e^{\mathrm i (m \varphi + k_z z)}

where :math:`k_z \in \mathbb R` and :math:`m \in \mathbb Z` are the parameters of the
solution (:func:`acoustotreams.scw_rPsi`, and :func:`acoustotreams.scw_Psi`). The radial part 
of the wave vector is defined as :math:`k_\rho = \sqrt{k^2 - k_z^2}` with the imaginary part 
of the square root to be taken non-negative. Note that here :math:`\rho` is a radial distance 
of cylindrical coordinates. The functions :math:`Z_m^{(n)}` are the Bessel and Hankel functions. 
For a complete set of solutions, it is necessary to select two of them. We generally use the (regular)
Bessel functions :math:`J_m = Z_m^{(1)}` and the Hankel functions of the first kind
:math:`H_m^{(1)} = Z_m^{(3)}`, which are singular at :math:`\rho \to 0` and correspond to radiating waves,
(:func:`acoustotreams.jv`, :func:`acoustotreams.hankel1`). The vector cylindrical waves are then

.. math::

    \mathbf{L}_{k_z, m}^{(n)}(k, \mathbf{r}) 
    = 
    \left[\frac{k_{\rho}}{k} Z'_{m}(k_{\rho}\rho) \hat{\boldsymbol{\rho}} 
    + \mathrm{i}\frac{m k_{\rho}}{k}\frac{Z_{m}(k_{\rho}\rho)}{k_{\rho}\rho} \hat{\boldsymbol{\varphi}} 
    + \frac{\mathrm{i} k_z}{k}Z_{m}(k_{\rho}\rho) \hat{\mathbf{z}} \right] 
    \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}

where we, again, normalized the functions by by :math:`k` in the medium (:func:`acoustotreams.vcw_rL`,
and :func:`acoustotreams.vcw_L`).

Spherical waves
---------------

Finally, we define the spherical wave solutions 

.. math::

    \psi_{lm}^{(n)}(k, \mathbf{r}) = z_l^{(n)}(kr) Y_{lm}(\theta, \phi)

given by (:func:`acoustotreams.ssw_rPsi`, and :func:`acoustotreams.ssw_Psi`)
where :math:`z_l^{(n)}` are the spherical Bessel and Hankel functions. 
We choose :math:`j_l = z_l^{(1)}` and :math:`h_l^{(1)} = z_l^{(3)}` in complete analogy 
to the cylindrical waves case (:func:`acoustotreams.spherical_jn`,
:func:`treams.special.spherical_hankel1`). :math:`Y_{lm}` are the spherical harmonics 
(:func:`treams.special.sph_harm`). The value :math:`l \in \mathbb N \cup \{0\}` refers to 
the angular momentum or degree. The projection of the angular momentum onto the z axis or order is 
:math:`m \in \mathbb Z` with :math:`|m| \leq l`. Hence, the vector spherical waves are defined as

.. math::

   \mathbf{L}_{lm}^{(n)} (k, \mathbf{r})
    = -\mathrm{i} \left[{z_l^{(n)}}' (kr) \mathbf{Z}_{lm}(\theta, \varphi) + \frac{z_l^{(n)}(kr)}{kr} \mathbf{Y}_{lm}(\theta, \varphi) \right]

(:func:`acoustotreams.vsw_rL`, and :func:`acoustotreams.vsw_L`) where

.. math::

    \mathbf X_{lm} (\theta, \varphi)
    = \mathrm i \sqrt{\frac{2 l + 1}{4 \pi l (l + 1)} \frac{(l - m)!}{(l + m)!}}
    \left(\mathrm i \pi_l^m(\cos\theta) \boldsymbol{\hat\theta}
    - \tau_l^m (\cos\theta) \boldsymbol{\hat\varphi}\right)
    \mathrm e^{\mathrm i m \varphi}
    \\
    \mathbf Y_{lm} (\theta, \varphi)
    = \mathrm i \sqrt{\frac{2 l + 1}{4 \pi l (l + 1)} \frac{(l - m)!}{(l + m)!}}
    \left(\tau_l^m (\cos\theta) \boldsymbol{\hat\theta}
    + \mathrm i \pi_l^m (\cos\theta) \boldsymbol{\hat\varphi}\right)
    \mathrm e^{\mathrm i m \varphi}
    \\
    \mathbf Z_{lm} (\theta, \varphi)
    = \mathrm i Y_{lm}(\theta, \varphi) \mathbf{\hat r}

are the vector spherical harmonics (:func:`acoustotreams.vsh_X`,
:func:`acoustotreams.vsh_Y`, and :func:`acoustotreams.vsh_Z`) imported from `treams.special`. 
These are themselves defined by the functions :math:`\pi_l^m(x) = \frac{m P_l^m(x)}{\sqrt{1 - x^2}}`,
:math:`\tau_l^m(x) = \frac{\mathrm d}{\mathrm d \theta}P_l^m(x = \cos\theta)`, and
the associated Legendre polynomials :math:`P_l^m` (:func:`acoustotreams.pi_fun`,
:func:`acoustotreams.tau_fun`, and :func:`acoustotreams.lpmv` are also imported from `treams.special`). 
The vector spherical harmonics are orthogonal to each other and normalized to 1 upon integration over the
solid angle.

Solutions to the acoustic master equations
==========================================

Up to now, we set up the acoustic master equations and found solutions to the scalar Helmholtz equation. 
The solutions to the acoustic master equations are then

.. math::

    p(k, \mathbf{r}) = \sum_\nu p_{\nu} \psi_\nu(k, \mathbf{r}) \\
    \mathbf{v}(k, \mathbf{r}) = \frac{-\mathrm{i}}{\rho c} \sum_\nu p_{\nu} \mathbf{L}_\nu(k, \mathbf{r})

where :math:`\nu` is just a placeholder for the actual parameters that index a concrete set of solutions.
The inverse prefactor :math:`Z = \rho c` is called acoustic impedance. 


References
==========

.. [1] H. Bruus, Acoustofluidics 2: Perturbation theory and 
   ultrasound resonance modes  Lab Chip 12, 20 (2012).
.. [2] L. D. Landau and E. M. Lifshitz, Theory of Elasticity: 
   Course of Theoretical Physics, Volume 7, 
   edited by E. M. Lifshitz, A. M. Kosevich, and  L. P. Pitaevskii
   (Butterworth-Heinemann, Oxford, England, UK, 1986).
.. [3] P. M. Morse and H. Feshbach, Methods of Theoretical Physics
   (McGraw-Hill, New York, 1953).