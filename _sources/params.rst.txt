.. testsetup::

   import numpy as np
   import acoustotreams

===========================================
Scalar basis sets and other core parameters
===========================================

Throughout the high-level functions and classes of *acoustotreams* a set of parameters appear
that define important underlying quantities for the calculation. First, these are the
different scalar basis sets that are used to solve scattering processes: the spherical,
cylindrical, and plane wave solutions. The word "scalar" emphasizes that the basis sets
do not have polarization. Closely related to these basis sets are the
mode types. The other parameters are the wave numbers in air and the materials as well as, 
in the case of calculations with periodicity involved, the lattice definitions and 
the phase shift between lattice sites.

Scalar basis sets
=================

As described in :doc:`theory` it is possible to solve linear-acoustics equations in different
coordinate systems. While being in principle equivalent, for different scenarios it is
beneficial to use suitable solution sets that represent the waves with sufficient
precision when truncated to a finite number of modes. The chosen finite number of
modes is given in the classes :class:`~acoustotreams.ScalarSphericalWaveBasis`,
:class:`~acoustotreams.ScalarCylindricalWaveBasis`, and :class:`~acoustotreams._core.ScalarPlaneWaveBasis`, which
are all children of the base call :class:`~acoustotreams._core.ScalarBasisSet`.

The modes of the spherical basis can are defined by their degree ``l``, and the order ``m``. 
The basis is then simply the collection of multiple of these modes, 
each given in a tuple with exactly that order, for example

.. doctest::

    >>> acoustotreams.ScalarSphericalWaveBasis([(0, 0), (1, -1), (1, 0), (1, 1)])
    ScalarSphericalWaveBasis(
        pidx=[0 0 0],
        l=[0 1 1 1],
        m=[0 -1  0  1],
        positions=[[0. 0. 0.]],
    )

results in a basis with four modes. The degree ``l`` is either 0 or 1, but the
order ``m`` is 0 for ``l=0`` and goes from -1 to 1 for ``l=1``. 
We see that there are also the fields ``pidx`` and ``positions``. This is 
a special case for the spherical (and later also the cyclindrical) wave basis. 
Sometimes, the fields are not expanded with respect to a
single point, but multiple positions. Then ``positions`` contains the their Cartesian
coordinates and ``pidx`` maps each mode to one of those coordinates. Here, the default
value of the expansion about a single origin is used. These basis sets behave mostly
like regular Python sets, we can check for example if a mode is in our basis set by

.. doctest::

    >>> (0, 1, 0) in acoustotreams.ScalarSphericalWaveBasis([(0, 0), (1, -1), (1, 0), (1, 1)])
    True

Equally, it is possible to use the regular comparisons and binary operators of Python
sets

.. doctest::

    >>> acoustotreams.ScalarSphericalWaveBasis([(0, 0), (1, -1), (1, 0), (1, 1)]) > {(0, 1, 0)}
    True
    >>> acoustotreams.ScalarSphericalWaveBasis([(0, 0), (1, -1), (1, 0), (1, 1)]) & {(0, 1, 0)}
    ScalarSphericalWaveBasis(
        pidx=[0],
        l=[1],
        m=[0],
        positions=[[0. 0. 0.]],
    )

However, because we want to use those basis sets later to index the rows and columns of
matrices, the order of the entries is fixed. Therefore, the equality operator is
stricter. Two basis sets are only considered equal when they have the same number modes
in the same order and the same positions.

.. doctest::

    >>> acoustotreams.ScalarSphericalWaveBasis([(0, 0), (1, 0)]) == acoustotreams.ScalarSphericalWaveBasis([(1, 0), (0, 0)])
    False

For convenience it is possible to create a default order up to a maximal multipolar
order

.. doctest::

    >>> acoustotreams.ScalarSphericalWaveBasis.default(2)
    ScalarSphericalWaveBasis(
        pidx=[0 0 0 0 0 0 0 0 0],
        l=[0 1 1 1 2 2 2 2 2],
        m=[ 0 -1  0  1 -2 -1  0  1  2],
        positions=[[0. 0. 0.]],
    )

where we now have a spherical wave basis up to the quadrupolar order.

The cyclindrical wave basis is mostly similar to the spherical basis. Instead of the
multipole degree ``l``, the z component of the wave vector ``kz`` is used

.. doctest::

    >>> acoustotreams.ScalarCylindricalWaveBasis([(.1, -1), (.1, 0), (.1, 1)])
    ScalarCylindricalWaveBasis(
        pidx=[0 0 0],
        kz=[0.1 0.1 0.1],
        m=[-1  0  1],
        positions=[[0. 0. 0.]],
    )

which is a real number. The default function takes a list of ``kz`` values and a maximal
absolute value of ``m``.

.. doctest::

    >>> acoustotreams.ScalarCylindricalWaveBasis.default([-.5, .5], 1)
    ScalarCylindricalWaveBasis(
        pidx=[0 0 0 0 0 0],
        kz=[-0.5 -0.5 -0.5  0.5  0.5  0.5],
        m=[-1  0  1 -1  0  1],
        positions=[[0. 0. 0.]],
    )

The cylindrical wave basis is particularly useful for systems with periodicity in the
z direction. Then, a basis with the diffraction orders up to a certain threshold can be obtained
by running

.. doctest::

    >>> acoustotreams.ScalarCylindricalWaveBasis.diffr_orders(kz=.1, mmax=1, lattice=2 * np.pi, bmax=1.05)
    ScalarCylindricalWaveBasis(
        pidx=[0 0 0 0 0 0 0 0 0],
        kz=[-0.9 -0.9 -0.9  0.1  0.1  0.1  1.1  1.1  1.1],
        m=[-1  0  1 -1  0  1 -1  0  1],
        positions=[[0. 0. 0.]],
    )

where ``bmax`` defines a distance in reciprocal space.

The plane wave basis behaves a little bit different. First, it is currently only defined
with respect to a single origin, so that the ``pidx`` and ``positions`` are not defined. Moreover,
the basis can be defined in two ways: :class:`ScalarPlaneWaveBasisByUnitVector` and
:class:`ScalarPlaneWaveBasisByComp`. In the first case, the definition is given by the unit
vector which, multiplied by the wave number in the medium, provides the full wave vector.
In the second case, two components of the wave vector are given and the remaining third
Cartesian component is defined such that it fulfils the dispersion relation.

. doctest::

    >>> acoustotreams.ScalarPlaneWaveBasisByUnitVector([(4, 0, 3)])
    ScalarPlaneWaveBasisByUnitVector(
        qx=[0.8],
        qy=[0.],
        qz=[0.6],
    )
    >>> acoustotreams.ScalarPlaneWaveBasisByComp([(1, 0)])
    ScalarPlaneWaveBasisByComp(
        kx=[1],
        ky=[0],
    )

By default, it is assumed, that the x and y components are given for the latter class,
but other components can also be chosen.

It is possible to convert between those basis sets by using the corresponding
functions

.. doctest::

    >>> pwbc = acoustotreams.ScalarPlaneWaveBasisByComp([(3, 0)])
    >>> pwbc.byunitvector(5)
    ScalarPlaneWaveBasisByUnitVector(
        qx=[0.6+0.j],
        qy=[0.+0.j],
        qz=[0.8+0.j],
    )
    >>> pwbuv = acoustotreams.ScalarPlaneWaveBasisByUnitVector([(0, 0, 1)])    
    >>> pwbuv.bycomp(1)
    ScalarPlaneWaveBasisByComp(
        kx=[0.],
        ky=[0.],
    )

Additionally, similar to the case of cylindrical waves, the basis by components can be
used for a range of diffraction orders

>>> acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders([0, 0], np.eye(2), 7)
    ScalarPlaneWaveBasisByComp(
        kx=[ 0.          0.          0.          6.28318531 -6.28318531],
        ky=[ 0.          6.28318531 -6.28318531  0.          0.        ],
    )


Polarization
============

In *acoustotreams* we consider acoustic waves that can be described by scalar pressure fields. 
The scalar waves are defined in :func:`acoustotreams.ssw_Psi`, :func:`acoustotreams.ssw_rPsi`,
:func:`acoustotreams.scw_Psi`, :func:`acoustotreams.scw_rPsi`, and :func:`acoustotreams.spw_Psi`.
These scalar waves correspond to the longitudinal vector fields: :func:`acoustotreams.vsw_L`, 
:func:`acoustotreams.vsw_rL`, :func:`acoustotreams.vcw_L`, :func:`acoustotreams.vcw_rL`, 
and :func:`acoustotreams.vpw_L`.

Mode types
==========

For some basis sets, there exist two different types of modes that distinguish
propagation features. For the spherical and cylindrical basis theses are `regular`
and `singular` modes. The former come through the use of (spherical) Bessel Functions
and the latter through the use of (spherical) Hankel functions of the first kind. The
regular modes are finite in the whole space. Thus, they are suitable for describing
incident modes or to expand a plane wave. The singular modes fulfil the radiation
condition and as such are used for the scattered fields.

For the plane wave basis of type (:class:`~acoustotreams.ScalarPlaneWaveBasisByComp`) only two
components of the wave vector are given and the third component is only implicitly
defined by the wave number and the material parameters. The application for this basis
is mostly within stratified media that are uniform or periodic in the two other
dimensions. Thus, the two given components of the wave vectors are conserved up to
reciprocal lattice vectors. To lift the ambiguity of the definition of the third
component, the mode types `up` and `down` are possible. They define, if the modes
propagate -- or decay for evanescent modes -- along the positive or negative direction
with respect to the third axis.

Wave number in air
==================

All calculations are executed in frequency domain. Instead of defining the frequency
:math:`\nu` or the angular frequency :math:`\omega` itself, *acoustotreams* works by using the
wave number in air

.. math::

    k_0 = \frac{2 \pi \nu}{c} = \frac{\omega}{c}

where :math:`c` is the speed of sound in air. In the code this real-valued number is
usually referred to by ``k0``. Implicitly, it is assumed throughout that all quantities,
like wave numbers, wave vectors, distances, or lattice vectors are given in the same
unit of (inverse) length.

Materials
=========

For materials there exists the class :class:`~acoustotreams.AcousticMaterial`, which holds the values
of the mass density, longitudinal speed of sound, and transverse speed of sound. The
default material is air (room temperature) and can be initialized without any parameters. For other cases,
the parameters can be given in the order above.

.. doctest::

    >>> acoustotreams.AcousticMaterial()
    AcousticMaterial(1.3, 343.0, 0.0)
    >>> acoustotreams.AcousticMaterial(1300, 2000, 1000)
    AcousticMaterial(1300, 2000, 1000)

It is also possible to calculate the impedance for longitudinal and transverse waves

.. doctest::

    >>> acoustotreams.AcousticMaterial().impedance
    '445.90000000000003'
    >>> acoustotreams.AcousticMaterial().impedancet
    '0.0'

Moreover, we can get the speeds from the Lamé parameters as

.. doctest::

    >>> acoustotreams.AcousticMaterial.from_params(1000, (1.6e7, 2e6))
    AcousticMaterial(1000, 141.421+0j, 44.7213+0j)

where the first is agrument is density and the second one is tuple :math:`(\lambda, \mu)`.

We can also get the speed of transverse waves from Poisson's ratio as

>>> acoustotreams.AcousticMaterial.from_pratio(1000, 141.421+0j, 0.44444)
    AcousticMaterial(1000, 141.421+0j, 44.7228+0j)


Lattices
========

The periodicity of arrangements is given by defining an instance of the class
:class:`~acoustotreams.Lattice`. A lattice can be one-, two-, or three-dimensional.

.. doctest::

    >>> acoustotreams.Lattice(1)
    Lattice(1.0, alignment='z')
    >>> acoustotreams.Lattice([[1, .5], [-.5, 1]])
    Lattice([[ 1.   0.5]
             [-0.5  1. ]], alignment='xy')
    >>> acoustotreams.Lattice([1, 2, 3])
    Lattice([[1. 0. 0.]
             [0. 2. 0.]
             [0. 0. 3.]], alignment='xyz')

The one- and two-dimensional lattices have to be aligned with one and two, respectively,
Cartesian axes. The default alignments are along the z axis for one-dimensional and in
the xy plane for the two-dimensional lattices. In the last example we see that it is
sufficient to just specify the diagonal entries. It is also possible to automatically
create lattices with special unit cells, for example

.. doctest::

    >>> acoustotreams.Lattice.hexagonal(2)
    Lattice([[2.    0.   ]
             [1.    1.732]], alignment='xy')

creates a hexagonal lattice with sidelength 2. It is also possible to extract a
lower-dimensional sublattice

.. doctest::

    >>> lat_3d = acoustotreams.Lattice([1, 2, 3])
    >>> acoustotreams.Lattice(lat_3d, "zx")
    Lattice([[0. 1.]
             [3. 0.]], alignment='zx')

or to combine and compare lattices

.. doctest::

    >>> acoustotreams.Lattice(1, "x") | acoustotreams.Lattice(2, "y")
    Lattice([[1. 0.]
             [0. 2.]], alignment='xy')
    >>> acoustotreams.Lattice([1, 2], "xy") & acoustotreams.Lattice([2, 3], "yz")
    Lattice(2.0, alignment='y')
    >>> acoustotreams.Lattice(1, "x") <= acoustotreams.Lattice([1, 2], "xy")
    True

The volume of the lattice can also be obtained

.. doctest::

    >>> acoustotreams.Lattice([[1, 0], [0, 1]]).volume
    1.0
    >>> acoustotreams.Lattice([[0, 1], [1, 0]]).volume
    -1.0

as we see the volume is "signed", i.e. it shows if the lattice vectors are in a
right-handed order, and the reciprocal lattice vectors can be computed

.. doctest::

    >>> acoustotreams.Lattice([1, 1]).reciprocal
    array([[ 6.283, -0.   ],
           [-0.   ,  6.283]])

Phase vector
============

The wave vector, often referred to as ``kpar``, specifies the phase relationship of
different lattice sites :math:`\mathrm{e}^{\mathrm i \mathbf{k}_\parallel \mathbf R}`.

.. doctest::

    >>> acoustotreams.WaveVector()
    WaveVector(nan, nan, nan)
    >>> acoustotreams.WaveVector(1)
    WaveVector(nan, nan, 1)
    >>> acoustotreams.WaveVector(1, "x")
    WaveVector(1, nan, nan)
    >>> acoustotreams.WaveVector((1, 2))
    WaveVector(1, 2, nan)
    >>> acoustotreams.WaveVector((1, 2, 3))
    WaveVector(1, 2, 3)

where unspecified directions are represented as ``nan``. The wave vectors can be
combined and compared.

.. doctest::

    >>> acoustotreams.WaveVector((1, 2)) | acoustotreams.WaveVector((2, 3), "yz")
    WaveVector(nan, 2, nan)
    >>> acoustotreams.WaveVector(1, "x") & acoustotreams.WaveVector(2, "y")
    WaveVector(1, 2, nan)
    >>> acoustotreams.WaveVector(1, "x") >= acoustotreams.WaveVector((1, 2))
    True

Note that the ordering is from less strict wave vector to the stricter one.