.. highlight:: python

.. only:: builder_html

===============================
Acoustic cylindrical T-Matrices
===============================

.. contents:: Table of contents
   :local:

Here, we cover acoustical cylindrical T-matrices, which are distinguished from the more
conventional (spherical) T-matrices through the use of scalar cylindrical waves as the
basis instead of scalar spherical waves. These waves are parametrized by the z component
of the wave vector :math:`k_z`, which describes their behavior in the z direction
:math:`\mathrm e^{\mathrm i k_z z}`, and the azimuthal order :math:`m`, which is also
used in scalar spherical waves.

The cylindrical T-matrices are suited for structures that are periodic in one dimension
(conventionally set along the z axis). Similarly to T-matrices of spheres that contain
the analytically known Mie coefficients, the cylindrical T-matrices for infinitely long
cylinders can also be calculated analytically.

Another similarity to spherical T-matrices are the possibilities to describe clusters of
objects in a local and global basis and to place these objects in a lattice. The
lattices can only extend in one and two dimensions; the z direction is implicitly
periodic already.

Infinitely long cylinders
=========================

The first simple object, for which we calculate the cylindrical T-matrix is an
infinitely long cylinder. Due to the rotation symmetry about the z axis this matrix 
is diagonal with respect to :code:`m` and due to the translation symmetry 
it is also diagonal with respect to :code:`kz`.


.. literalinclude:: examples/cylinder_tmatrixc.py
    :language: python
    :lines: 7-13

For such infinitely long structures it makes more sense to talk about cross width
instead of cross section. We obtain the averaged scattering and extinction cross width
by 

.. literalinclude:: examples/cylinder_tmatrixc.py
    :language: python
    :lines: 15-16

Again, we can also select specific modes only, for example the modes with :math:`m = 0` and :math:`m = \pm 1`.

.. literalinclude:: examples/cylinder_tmatrixc.py
    :language: python
    :lines: 18-21

to calculate their cross width. Evaluating the field intensity in the xy plane is similar to the case of a sphere.

.. literalinclude:: examples/cylinder_tmatrixc.py
    :language: python
    :lines: 23-42

For the infinite cylinder we can also calculate the radiation pattern as a function of the angle :math:`\varphi`.
The red arrow indicates the direction of incidence of the plane wave.

.. plot:: examples/cylinder_tmatrixc.py

Cylindrical T-matrices for one-dimensional arrays of spherical T-matrices
=========================================================================

For the next example let us look at the system of spheres on a one-dimensional
lattice again (:ref:`tmatrix:One-dimensional arrays (along z)`). They fulfil all
properties that define structures where the use of cylindrical waves is beneficial,
namely they have a finite extent in the xy plane and they are periodic along the
z direction. 

The initial setup of our calculation begins with spheres in the spherical wave
basis and placing them in a chain. This is the same procedure as in
:ref:`tmatrix:One-dimensional arrays (along z)`.

.. literalinclude:: examples/chain_tmatrixc.py
    :language: python
    :lines: 

Next, we convert this chain in the spherical wave basis to a suitable cylindrical wave
basis.

.. literalinclude:: examples/chain_tmatrixc.py
    :language: python
    :lines: 

We chose to add the first three diffraction orders (plus a 0.1 margin to avoid problems
with floating point comparisons).

Further, we set up the illumination and calculate the scattering by the usual
procedure.

.. literalinclude:: examples/chain_tmatrixc.py
    :language: python
    :lines: 

We evaluate the fields in two regions. Outside of the circumscribing cylinders we can
use the fast cylindrical wave expansion. Inside of the circumscribing cylinders but
outside of the spheres we can use the method of
:ref:`tmatrix:One-dimensional arrays (along z)`.

Finally, we can plot the results. To illustrate the periodicity better, three unit cells
are shown.

.. plot:: examples/chain_tmatrixc.py

One-dimensional arrays (along the x axis)
=========================================

Two-dimensional arrays (in the xy plane)
=========================================