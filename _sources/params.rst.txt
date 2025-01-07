.. testsetup::

   import numpy as np
   import acoutreams

===========================================
Scalar basis sets and other core parameters
===========================================

Throughout the high-level functions and classes of *acoutreams* a set of parameters appear
that define important underlying quantities for the calculation. First, these are the
different scalar basis sets that are used to solve scattering processes: the spherical,
cylindrical, and plane wave solutions. The word "scalar" emphasizes that the basis sets
do not have polarization. Closely related to these basis sets are the
mode types. The other parameters are the vacuum wave numbers and the materials as well as, 
in the case of calculations with periodicity involved, the lattice definitions and 
the phase shift between lattice sites.

Scalar basis sets
=================

As described in :doc:`theory` it is possible to solve linear-acoustics equations in different
coordinate systems. While being in principle equivalent, for different scenarios it is
beneficial to use suitable solution sets that represent the waves with sufficient
precision when truncated to a finite number of modes. The chosen finite number of
modes is given in the classes :class:`~acoutreams.ScalarSphericalWaveBasis`,
:class:`~acoutreams.ScalarCylindricalWaveBasis`, and :class:`~acoutreams._core.ScalarPlaneWaveBasis`, which
are all children of the base call :class:`~acoutreams._core.ScalarBasisSet`.

The modes of the spherical basis can are defined by their degree ``l``, and the order ``m``. 
The basis is then simply the collection of multiple of these modes, 
each given in a tuple with exactly that order, for example

.. doctest::

    >>> acoutreams.ScalarSphericalWaveBasis([(0, 0), (1, -1), (1, 0), (1, 1)])
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

    >>> (0, 1, 0) in acoutreams.ScalarSphericalWaveBasis([(0, 0), (1, -1), (1, 0), (1, 1)])
    True

Equally, it is possible to use the regular comparisons and binary operators of Python
sets

.. doctest::

    >>> acoutreams.ScalarSphericalWaveBasis([(0, 0), (1, -1), (1, 0), (1, 1)]) > {(0, 1, 0)}
    True
    >>> acoutreams.ScalarSphericalWaveBasis([(0, 0), (1, -1), (1, 0), (1, 1)]) & {(0, 1, 0)}
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

    >>> acoutreams.ScalarSphericalWaveBasis([(0, 0), (1, 0)]) == acoutreams.ScalarSphericalWaveBasis([(1, 0), (0, 0)])
    False

For convenience it is possible to create a default order up to a maximal multipolar
order

.. doctest::

    >>> acoutreams.ScalarSphericalWaveBasis.default(2)

where we now have a spherical wave basis up do quadrupolar order.