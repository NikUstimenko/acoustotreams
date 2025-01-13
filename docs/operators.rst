=========
Operators
=========

There are numerous operators implemented in *acoustotreams*. They replicate to some extend the
way active transformations with operators work on linear mapping :math:`M`

.. math::

    O(x) M O^{-1}(x)

where :math:`O` is the transformation we want to apply and and :math:`x` is a
parameter of that transformation. Similarly,

.. math::

    O(x) \psi

transforms state :math:`\psi` that can be represented by a vector. We attempt to
replicate this transformation notation in code and to extend it by other useful
functions. Later, the state :math:`\psi` can often be the expansion coefficients of
a wave, the linear mapping could be the T-matrix and the transformation operator could
be a rotation.

Each operator is implemented as a class. The class is instantiated with the parameter
:code:`op = Operator(x)`. At this stage it is only an abstract operator. A concrete
representation can be obtained by calling the operator with the necessary keyword
arguments, e.g. :code:`op(basis=concrete_basis)`, which will return a array-like
structure.

However, to be able to replicate the above mathematical notation, it is also possible to
use the matrix multiplication operator between an array and the operator. The array
needs to have the necessary keywords as attributes, this is for what :doc:`acousticsarray`
come in handy. For such an array :code:`arr` it is possible to type :code:`op @ arr` or
:code:`op @ arr @ op.inv`. The inverse of operators is implemented for many operators
but some have no inverse defined.

Sometimes, it can come in handy to directly apply an operator to an array without
defining the abstract operator before. This can be achieved by
:code:`arr.operator.apply_left(x)` and :code:`arr.operator.apply_right(x)`, which are
equivalent to :code:`op @ arr` and :code:`arr @ op.inv`, respectively. The function
:code:`arr.op()` is also defined. For arrays with ``ndim`` equal to 1 or without an
inverse, it is equivalent to :code:`arr.op.apply_left()`, otherwise it corresponds to
:code:`op @ arr @ op.inv`.

Rotation
========

As already mentioned, a common transformation are rotations. The representation of the
rotation operator depends on which basis we use. For example, for the T-matrix using the
spherical wave basis, such a rotation is represented by the Wigner D-matrix elements,
but for plane waves this would look different. In *acoustotreams* we can create such an
abstract rotation operator by using :class:`~acoustotreams.Rotate`

.. doctest::

    >>> r = acoustotreams.Rotate(np.pi)

which can then be converted to a representation by calling it with the basis argument.

.. doctest::

    >>> r(basis=acoustotreams.ScalarSphericalWaveBasis.default(1))
    AcousticsArray(
        [[...]],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0 -1 0 1],
        positions=[[0. 0. 0.]],
    ),
    )

If it is multiplied with an array that defines the attribute `basis` it will
automatically take that attribute.

.. doctest::

    >>> t = acoustotreams.AcousticsArray(np.eye(4), basis=acoustotreams.ScalarSphericalWaveBasis.default(1))
    >>> r @ t @ r.inv
    AcousticsArray(
        [[...]],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0 -1 0 1],
        positions=[[0. 0. 0.]],
    ),
    )

Here, we also use the property `inv` to get the inverse rotation. Moreover, for
instances of :class:`~acoustotreams.AcousticsArray`, we can get the same result by calling the
correspondingly named attribute

    >>> phi = 1
    >>> r = acoustotreams.Rotate(phi)
    >>> (r @ t @ r.inv == t.rotate(phi)).all()
    True

which also has the methods `apply_left` and `apply_right` to only apply the operator
from one side. For some basis sets only rotations about the z axis are possible, while
other basis sets allow rotations including all three Euler angles.

Translation
===========

The next transformation that is implemented are translations where the parameter is the
Cartesian translation vector.

.. doctest::

    >>> t = acoustotreams.AcousticsArray(np.eye(6), basis=acoustotreams.ScalarSphericalWaveBasis.default(1), k0=1)
    >>> t.translate([1, 2, 3])
    AcousticsArray(
        [[...]],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0 -1 0 1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        material=Material(1.3, 343.0, 0.0),
    )

For the translation we have to specify the basis and the vacuum wave number. In the
result we can see that the default material of the embedding is air.

.. note::

    The rotation and translation operators applied to a spherical or cylindrical basis
    with multiple positions, will rotate or translate each position independently from
    the others. This results in block-diagonal matrices with respect to the different
    positions in such a case.

Expand in a different basis
===========================

The expansion in a different basis set is a little bit more complicated due to the
number of possible combinations of which basis set can be expanded in which other basis
sets. Therefore, we will treat each source basis set separately in the following.

Also, here the notion of abstract operator and concrete representation breaks down to
some extent because it makes little sense to first define an abstract expansion in,
e.g., spherical waves without specifying the relevant multipoles. Thus, the concrete
representation of the target basis is the argument of the operator.

Plane waves
-----------

Plane waves can be expanded into a different set of plane waves and into regular
spherical and cylindrical waves. The expansion into a different set of plane waves
is basically just a matching of the wave vectors.

.. doctest::

    >>> plw = acoustotreams.plane_wave_scalar([0, 3, 4], k0=5, material=1.3)
    >>> acoustotreams.Expand(acoustotreams.ScalarPlaneWaveBasisByComp.default([[0, 3]])) @ plw
    AcousticsArray(
        [1.+0.j],
        basis=ScalarPlaneWaveBasisByComp(
        kx=[0.],
        ky=[3.],
    ),
        k0=5.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='up',
    )

For example, here we change from the expansion in
:class:`~acoustotreams.ScalarPlaneWaveBasisByUnitVector` to the expansion by x and y components.
For such a basis change, it is necessary that the material and the wave number is
specified.

Next, we can expand this plane wave also in cylindrical and in spherical waves.

.. doctest::

    >>> acoustotreams.Expand(acoustotreams.ScalarCylindricalWaveBasis.default([4], 1)) @ plw
    AcousticsArray(
        [1.+0.j 1.+0.j 1.+0.j],
        basis=ScalarCylindricalWaveBasis(
        pidx=[0 0 0],
        kz=[4. 4. 4.],
        m=[-1 0 1],
        positions=[[0. 0. 0.]],
    ),
        k0=5.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='regular',
    )
    >>> acoustotreams.Expand(acoustotreams.ScalarSphericalWaveBasis.default(1)) @ plw
    AcousticsArray(
        [ 3.5449077 +0.00000000e+00j -2.60496452+1.59508073e-16j
          0.        +4.91196820e+00j -2.60496452-1.59508073e-16j   ],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[1 1 1 1],
        m=[0 -1 0 1],
        positions=[[0. 0. 0.]],
    ),
        k0=5.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='regular',
    )

Spherical waves
---------------

Next, we expand spherical waves. In comparison to the plane waves, spherical waves have
the added difficulty of the categorization of "regular" and "singular" functions and the
distinction of global and local basis sets.

In a simple case, we want to expand a spherical wave that is centered not at the origin
and expand it around the origin

.. doctest::

    >>> off_centered_swb = acoustotreams.ScalarSphericalWaveBasis.default(1, positions=[[1, 0, 0]])
    >>> sssw = acoustotreams.spherical_wave_scalar(1, 0, basis=off_centered_swb, k0=1, material=1.3, modetype="singular")
    >>> ex = acoustotreams.Expand(acoustotreams.ScalarSphericalWaveBasis.default(1))
    >>> ex @ sssw
    AcousticsArray(
        [ 0.+0.j 0.+0.j 0.90350604+0.j 0.+0.j ],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0 -1 0 1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='singular',
    )

We defined the wave as a singular wave and, if nothing is explicitly specified, the
expansion into other spherical waves is taken as the same type of field. So, a singular
field will be expanded again in singular modes and a regular field is expanded in
regular modes. However, we can also change the type of mode, when the field is expanded
around a different origin

.. doctest::

    >>> ex = acoustotreams.Expand(acoustotreams.ScalarSphericalWaveBasis.default(1), "regular")
    >>> ex @ sssw
    AcousticsArray(
        [ 0.+0.j 0.+0.j 0.90350604-4.14531987j 0.+0.j ],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0 -1 0 1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='regular',
    )

for this we had to define the ``modetype`` for the expand operator.

Next, we want to look at the expansion of a global field into a local field at multiple
origins, which works quite similarly

.. doctest::

    >>> ssw_global = acoustotreams.spherical_wave_scalar(1, 0, k0=1, material=1.3, modetype="regular")
    >>> local_swb = acoustotreams.ScalarSphericalWaveBasis.default(1, 2, positions=[[0, 0, 1], [0, 0, -1]])
    >>> ssw_global.expand.apply_left(local_swb)
    AcousticsArray(
        [ 0.52163945+0.j  0.+0.j  0.71740088+0.j  0.+0.j
         -0.52163945+0.j  0.+0.j  0.71740088+0.j  0.+0.j],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0 1 1 1 1],
        l=[0 1 1 1 0 1 1 1],
        m=[0 -1 0 1 0 -1 0 1],
        positions=[[ 0.  0.  1.], [ 0.  0. -1.]],
    ),
        k0=1.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='regular',
    )

For the translations within only regular or only singular waves it is possible to
expand back into the same basis set in this case corresponds to the multiplication by a
unit matrix.

.. doctest::

    >>> sw_global.expand(acoustotreams.ScalarSphericalWaveBasis.default(1))
    AcousticsArray(
        [0.+0.j 0.+0.j 1.+0.j 0.+0.j],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0 -1 0 1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='regular',
    )

For translations from singular to regular waves, the same basis set means that a
zero matrix is returned.

.. doctest::

    >>> ssw_global.expand.apply_right(acoustotreams.ScalarSphericalWaveBasis.default(1), "singular")
    AcousticsArray(
        [0.+0.j 0.+0.j 0.+0.j 0.+0.j],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0 -1 0 1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='singular',
    )

Besides that the expansion of spherical waves in different basis sets results in dense
matrices.

The expansion of spherical waves into cylindrical or plane waves is a continuous
spectrum and is currently not implemented.

Cylindrical waves
-----------------

Cylindrical waves are similar to spherical waves, in the sense, that they can be
separated into regular and singular modes and that they can be defined with multiple
origins within *acoustotreams*. Therefore, the expansion within cylindrical waves follows the
same properties than spherical waves.

.. doctest::

    >>> off_centered_cwb = acoustotreams.ScalarCylindricalWaveBasis.default([0], 1, positions=[[1, 0, 0]])
    >>> sscw = acoustotreams.cylindrical_wave_scalar(0, 1, basis=off_centered_cwb, k0=1, material=1.3, modetype="singular")
    >>> ex =  acoustotreams.Expand(acoustotreams.ScalarCylindricalWaveBasis.default([0], 1))
    >>> ex @ sscw
    AcousticsArray(
        [0.11490348-1.40716185e-17j -0.44005059+5.38906541e-17j 0.76519769+0.00000000e+00j],
        basis=ScalarCylindricalWaveBasis(
        pidx=[0 0 0],
        kz=[0. 0. 0.],
        m=[-1  0  1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='singular',
    )

Additionally, it is possible to expand a cylindrical wave into spherical waves. Note,
that waves defined with multiple origins get each expanded separately. The positions
of the spherical and cylindrical waves must be equal.

.. doctest::

    >>> sscw = acoustotreams.cylindrical_wave_scalar(0, 1, k0=1, material=1, modetype="regular")
    >>> sscw.expand(acoustotreams.ScalarSphericalWaveBasis.default(1))
    AcousticsArray(
        [0.+0.j  0.+0.j  0.+0.j -4.34160753+0.j],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0  -1  0  1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='regular',
    )

The inverse of this expansion is not implemented.

The expansion of cylindrical waves into plane waves is a continuous spectrum and is not
implemented.

Expand in a different basis with periodic boundaries
====================================================

There is a special case of expansion implemented for the case of periodic boundaries
when using spherical or cylindrical waves. These expansions are needed to compute the
acoustic interaction between particles within a lattice. It is assumed that the
given basis with singular modes are repeated periodically in the given lattice
structure. Then, these fields are expanded as regular fields in a single unit cell.

.. doctest::

    >>> sscw = acoustotreams.cylindrical_wave_scalar(0, 1, k0=1, material=1.3, modetype="singular")
    >>> sscw.expandlattice(1, 0)
    AcousticsArray(
        [2.-3.8655259j  0.+0.j  1.+1.23397896j],
        basis=ScalarCylindricalWaveBasis(
        pidx=[0 0 0],
        kz=[0. 0. 0.],
        m=[-1  0  1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        kpar=WaveVector(0, nan, 0.0),
        lattice=Lattice(1.0, alignment='x'),
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='regular',
    )
    >>> sssw = acoustotreams.spherical_wave_scalar(1, 0, k0=1, material=1.3)
    >>> sssw.expandlattice([1, 2], [0, 0])
    AcousticsArray(
        [0.+0.j        0.+0.j        8.42477796-8.5237243j
         0.+0.j       ],
        basis=ScalarSphericalWaveBasis(
        pidx=[0 0 0 0],
        l=[0 1 1 1],
        m=[0 -1  0  1],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        kpar=WaveVector(0, 0, nan),
        lattice=Lattice([[1. 0.]
             [0. 2.]], alignment='xy'),
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='regular',
    )

The inverse of this operator is not implemented. Additionally, it is possible to expand
the periodic field into a different basis set. Spherical waves in a one-dimensional
lattice along the z axis can be expanded in cylindrical waves

.. doctest::

    >>> sssw = acoustotreams.spherical_wave_scalar(1, 0, k0=1, material=1.3)
    >>> ex = acoustotreams.ExpandLattice(basis=acoustotreams.ScalarCylindricalWaveBasis.diffr_orders([.1], 0, 7, 1))
    >>> ex @ sssw
    AcousticsArray(
        [0.+0.17490069j 0.-0.02192843j 0.-0.21875755j],
        basis=ScalarCylindricalWaveBasis(
        pidx=[0 0 0 0 0 0],
        kz=[-0.7975979  0.1        0.9975979],
        m=[0 0 0 0 0 0],
        positions=[[0. 0. 0.]],
    ),
        k0=1.0,
        kpar=WaveVector(nan, nan, 0.1),
        lattice=Lattice(7.0, alignment='z'),
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='singular',
    )

where the lattice and the wave vector are implicitly defined by the use of the
class method :func:`acoustotreams.ScalarCylindricalWaveBasis.diffr_orders`. Similarly, spherical
waves in a two-dimensional lattice in the xy plane can be expanded in plane waves.

.. doctest::

    >>> ex = acoustotreams.ExpandLattice(basis= acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders([.1, 0], [7, 7], 1))
    >>> ex @ sssw
    AcousticsArray(
        [0.-0.06265266j 0.-0.06265266j 0.-0.06265266j 0.-0.06265266j 0.-0.06265266j],
        basis=ScalarPlaneWaveBasisByComp(
        kx=[ 0.1        0.1        0.1        0.9975979 -0.7975979],
        ky=[ 0.         0.8975979 -0.8975979  0.         0.       ],
    ),
        k0=1.0,
        kpar=WaveVector(0.1, 0, nan),
        lattice=Lattice([[7. 0.]
             [0. 7.]], alignment='xy'),
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='up',
    )

Cylindrical waves, that themselves are periodic in the z direction, in a one-dimensional
lattice along the x axis can also be expanded in plane waves.

.. doctest::

    >>> sscw = acoustotreams.cylindrical_wave_scalar(0, 1, k0=1, material=1.3)
    >>> ex = acoustotreams.ExpandLattice(basis=acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders([0, .1], acoustotreams.Lattice([7, 7], "zx"), 1))
    >>> ex @ sscw
    AcousticsArray(
        [0.28571429-0.02871537j 0.28571429-4.1146983j  0.28571429+0.37780019j
         0.        +0.j         0.        +0.j        ]
        basis=ScalarPlaneWaveBasisByComp(
        kz=[ 0.         0.         0.         0.8975979 -0.8975979],
        kx=[ 0.1        0.9975979 -0.7975979  0.1        0.1      ],
    ),
        k0=1.0,
        kpar=WaveVector(nan, nan, 0.0),
        lattice=Lattice(7.0, alignment='x'),
        material=AcousticMaterial(1.3, 343.0, 0.0),
        modetype='up',
    )

Permute the axes
================

The permute operator is only implemented for plane waves.
For this type of waves, the rotation is only implemented about the z axis. These
rotations then do not include a relabeling of the Cartesian axes, for example
:math:`(x', y', z') = (z, x, y)`. This operation is implemented separately as
permutation, meaning the axes labels get permuted.

.. doctest::

    >>> plw = acoustotreams.plane_wave_scalar([2, 3, 6])
    >>> plw
    AcousticsArray(
        [1.+0.j],
        basis=ScalarPlaneWaveBasisByUnitVector(
        qx=[0.286],
        qy=[0.429],
        qz=[0.857],
    ),
    )
    >>> plw.permute()
    AcousticsArray(
        [1.+0.j],
        basis=ScalarPlaneWaveBasisByUnitVector(
        qx=[0.857],
        qy=[0.286],
        qz=[0.429],
    ),
    )

Evaluate the field
==================

From a programming perspective, the evaluation of the field values at specified points
is also implemented by a couple of operators. The pressure field :math:`p`, and
the velocity field :math:`\boldsymbol v` in Cartesian coordinates can be computed.

.. doctest::

    >>> sssw = acoustotreams.spherical_wave_scalar(1, 0, k0=1, material=1.3, modetype="regular")
    >>> sssw.pfield([0, 0, 1])
    AcousticsArray(
        [0.14715177+0.j],
    )
    >>> sssw.vfield([0, 0, 1])
    AcousticsArray(
        [0.+0.j         0.+0.j         0.-0.00026203j],
    )

Evaluate the far-field amplitude
================================

Similarly, we can also calculate the far-field amplitudes of the scattered pressure and velocity fields
defined as :math:`p = \frac{\mathrm e^{\mathrm i k r}}{r}p_0` and :math:`\boldsymbol v = \frac{\mathrm e^{\mathrm i k r}}{r}\boldsymbol v_0` for spherical waves,
and :math:`p = \frac{\mathrm e^{\mathrm i k_{\rho} \rho}}{\sqrt{\rho}}p_0` and :math:`\boldsymbol v = \frac{\mathrm e^{\mathrm i k_{\rho} \rho}}{\sqrt{\rho}}\boldsymbol v_0`
for cylindrical waves. For a singular spherical wave, the amplitudes :math:`p_0` and :math:`\boldsymbol v_0` are calculated
using :func:`acoustotreams.ssw_psi` and :func:`acoustotreams.ssw_l`, respectively; for a singular cylindrical wave, 
using :func:`acoustotreams.scw_psi` and :func:`acoustotreams.scw_l`.

.. doctest::

    >>> sssw = acoustotreams.spherical_wave_scalar(1, 0, k0=1, material=1.3, modetype="singular")
    >>> sssw.pamplitudeff([0, 0, 1])
    AcousticsArray(
        [-0.48860251+0.j],
    )
    >>> sssw.vamplitudeff([0, 0, 1])
    AcousticsArray(
        [-0.00109577+0.j  0.+0.j  0.+0.j],
    )

Note that :math:`v_0` is computed in spherical coordinates for the spherical wave basis and
cylindrical coordinates for the cylindrical wave basis.