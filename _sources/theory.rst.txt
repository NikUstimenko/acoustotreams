======
Theory
======

These sections give an overview over theory that underlies the different aspects of acoustotreams.
For more in-depth information, there is also a list of associated :ref:`about:Publications`.

Master equations of linear monochromatic acoustics
==================================================

The master equations that describe propagation of pressure waves can be written 
in frequency domain in the absence of external sources and shear (transverse) 
waves as [1]_

.. math::

    \begin{pmatrix}
    \boldsymbol{\nabla}& 0 \\
    0& \boldsymbol{\nabla} \cdot
    \end{pmatrix}
    \begin{pmatrix}
    p\\
    \boldsymbol{v}
    \end{pmatrix} = 
    \mathrm i \omega
    \begin{pmatrix}
    0&  \rho\\
    \frac{1}{\rho c^2}& 0
    \end{pmatrix}
    \begin{pmatrix}
    p\\
    \boldsymbol{v}
    \end{pmatrix}

where :math:`p = p(\boldsymbol{r}, \omega)` and :math:`\boldsymbol v = \boldsymbol v(\boldsymbol r, \omega)`
are the pressure and velocity fields (:func:`acoustotreams.pfield`, :func:`acoustotreams.vfield`).
All these quantities are complex valued fields, that depend on the angular frequency :math:`\omega` and 
the position :math:`\boldsymbol r`. The material parameters are the speed of sound :math:`c`, and
the mass density :math:`\rho`. Conventionally, within acoustotreams the (air) wave number
:math:`k_0 = \frac{\omega}{c_0}` is generally used to express the frequency. Here :math:`c_0` is
the speed of sound in air being 343 m/s.

In acoustotreams we can also find the pressure and velocity fields in the presence of objects
that support both pressure and shear waves while the background media support *no* shear waves. 
In this case we explicitly solve the equation  for the displacement field :math:`\boldsymbol u` [2]_

.. math::

    \omega^2 \boldsymbol u + c^2 \boldsymbol \Delta \boldsymbol u + (c^2 - c_t^2) \boldsymbol \nabla (\boldsymbol \nabla \cdot \boldsymbol u) = 0

where :math:`\rho` is the mass density, :math:`c` is the speed of longitudinal pressure waves, 
and :math:`c_t` is the speed of transverse shear waves (:func:`acoustotreams.AcousticMaterial`).
Due to the requirement of isotropy these quantities are all scalar. In a medium with :math:`c_t \equiv 0`
the last equation coincides with the second equation from the first set.

For the transformation to the time domain we use for a general function
:math:`f(\omega)`

.. math::

    f(t) = \int_{-\infty}^\infty \mathrm d t f(\omega) \mathrm e^{-\mathrm i \omega t}

as Fourier transformation convention, and thus the inverse transformation is

.. math::

    f(\omega)
    = \int_{-\infty}^\infty \frac{\mathrm d \omega}{2 \pi}
    f(t) \mathrm e^{\mathrm i \omega t}


References
==========

.. [1] H. Bruus, Acoustofluidics 2: Perturbation theory and 
   ultrasound resonance modes  Lab Chip 12, 20 (2012).
.. [2] L. D. Landau and E. M. Lifshitz, Fluid Mechanics: Landau
   and Lifshitz: Course of Theoretical Physics, Volume 6 (Elsevier,
   Amsterdam, 2013).