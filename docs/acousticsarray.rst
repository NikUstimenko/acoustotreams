=============================
Acoustics-aware arrays
=============================

.. testsetup::

   import numpy as np
   import acoustotreams

A core building block of the underlying features of *acoustotreams* are acoustics-aware arrays.
In most of their properties they behave similar to numpy arrays and one can easily
change the type and mix them

.. doctest::

    >>> np.array([1, 2]) * acoustotreams.AcousticsArray([2, 3])
    AcousticsArray(
        [2, 6],
    )
    >>> np.array([1, 2]) @ acoustotreams.AcousticsArray([2, 3])
    8

but they have mainly two features added. First, they derive from
:class:`treams.util.AnnotatedArray`, so that they can carry annotations with them, but these
annotations are restricted to the physical quantities (as described in :doc:`params`).
Secondly, they offer special methods to create matrices for common transformations like
rotations and translations, which are described in more detail in :doc:`operators`.

Special properties
==================

.. doctest::

    >>> acoustotreams.AcousticsArray([[0, 1], [2, 3]], k0=(1, 2))
    AcousticsArray(
        [[0, 1],
         [2, 3]],
        k0=(1.0, 2.0),
    )

In this example you can notice that the values for the air wave number ``k0`` were
converted from integers to floats. Thus, trying to use
:code:`acoustotreams.AcousticsArray([1], k0=1j)` will raise an error, because the complex number
cannot be interpreted as a float. Additional special keywords are `basis`, `kpar`,
`lattice`, `material`, and `modetype`. These properties can also be accessed
by setting the corresponding attribute

.. doctest:: 

    >>> m = acoustotreams.AcousticsArray([1, 2])
    >>> m.material = [4000, 1000]
    >>> m
    AcousticsArray(
        [1, 2],
        material=AcousticMaterial(4000, 1000, 0),
    )

where we now have a material with the mass density :math:`4000 kg/m^3` and speed of sound :math:`1000 m/s`. 
As with its parent class these properties are also compared and merged when using operations on these objects

.. doctest::

    >>> acoustotreams.AcousticsArray([0, 1], k0=1) + acoustotreams.AcousticsArray([2, 3], material=[2000, 1000, 0])
    AcousticsArray(
        [2, 4],
        k0=1.0,
        material=AcousticMaterial(2000, 1000, 0),
    )

and using conflicting values will raise a warning, for example
:code:`acoustotreams.AcousticsArray([0, 1], k0=1) + acoustotreams.AcousticsArray([2, 3], k0=2)`
emits :code:`treams/util.py:249: AnnotationWarning: at index 0: overwriting key 'k0'`.
The special properties have also a unique behavior when appearing in matrix
multiplications. If one of the two matrices has the special property not set, it becomes
"transparent" to it. Check out the difference between

.. doctest::

    >>> np.ones((2, 2)) @ acoustotreams.AcousticsArray([1, 2], k0=1.0)
    AcousticsArray(
        [3., 3.],
        k0=1.0,
    )

and 

.. doctest::

    >>> np.ones((2, 2)) @ treams.util.AnnotatedArray([1, 2], k0=(1.0,))
    AnnotatedArray(
        [3., 3.],
        AnnotationSequence(AnnotationDict({})),
    )

where besides the obvious difference in array types, the property `k0` is preserved.

The full list of special properties is:

======== ==============================================================
Name     Description
======== ==============================================================
basis    Basis set: spherical, cylindrical, planar
k0       Wave number in air (at the room temperature)
kpar     Phase relation in lattices (:class:`acoustotreams.WaveVector`)
lattice  Definition of a lattice (:class:`acoustotreams.Lattice`)
modetype Modetype, depends on wave (:ref:`params:Mode types`)
material Embedding material (:class:`acoustotreams.AcousticMaterial`)
======== ==============================================================