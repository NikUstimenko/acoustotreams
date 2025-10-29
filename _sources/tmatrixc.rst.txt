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

For our next example we want to look at the system of spheres on a one-dimensional
lattice again (:ref:`tmatrix:One-dimensional arrays (along z)`). They fulfil all
properties that define structures where the use of cylindrical waves is beneficial,
namely they have a finite extent in the xy plane and they are periodic along the
z direction.

Thus, the initial setup of our calculation starts with spheres in the spherical wave
basis and placing them in a chain. This is the same procedure as in
:ref:`tmatrix:One-dimensional arrays (along z)`. We define a nonzero component of 
the wave vector along the z axis :math:`k_z` in the air. 

.. literalinclude:: examples/chain_tmatrixc.py
    :language: python
    :lines: 7-19

Here we also define :math:`n` as a "refractive index" for the background medium of the spheres
with respect to the air and divide :math:`k_z` by :math:`n` when it is needed. Next, we convert 
this chain in the spherical wave basis to a suitable cylindrical wave basis.

.. literalinclude:: examples/chain_tmatrixc.py
    :language: python
    :lines: 21-23

We chose to add the first three diffraction orders (plus a 0.1 margin to avoid problems
with floating point comparisons).

Finally, we set-up the incident wave and calculate the scattering coeffcients by the usual
procedure.

.. literalinclude:: examples/chain_tmatrixc.py
    :language: python
    :lines: 25-31

We evaluate the fields in two regions. Outside of the circumscribing cylinders we can
use the fast cylindrical wave expansion. Inside of the circumscribing cylinders but
outside of the spheres we can use the method of
:ref:`tmatrix:One-dimensional arrays (along z)`.

.. literalinclude:: examples/chain_tmatrixc.py
    :language: python
    :lines: 33-52

Finally, we can plot the results. To illustrate the periodicity better, three unit cells
are shown.

.. plot:: examples/chain_tmatrixc.py


Clusters
========

Similarly to the case of spheres we can also calculate the response from a cluster of
objects. As an example let us simulate a cylinder together with a chain of spheres
in the cylindrical wave basis as described in the previous section.

Hence, we set up first the spheres in the chain and convert them to the cylindrical wave
basis as before

.. literalinclude:: examples/cluster_tmatrixc.py
    :language: python
    :lines: 7-22

Then, we create the T-matrix of the cylinder in the cylindrical wave basis.

.. literalinclude:: examples/cluster_tmatrixc.py
    :language: python
    :lines: 24

Finally, we construct the cluster, find the T-matrix of the interacting system,
and then scattered field coeffcients for the incident plane wave.

.. literalinclude:: examples/cluster_tmatrixc.py
    :language: python
    :lines: 26-32

The scattered pressure field within three unit cells is shown below

.. plot:: examples/cluster_tmatrixc.py

One-dimensional arrays (along the x axis)
=========================================

Now, we take the chain of spheres and cylinder and place them in a grating structure
along the x direction. We start again by defining the parameters and calculating the 
relevant cylindrical T-matrices.

.. literalinclude:: examples/grating_tmatrixc.py
    :language: python
    :lines: 7-25

Next, we create the cluster and, as usual, let it interact within a lattice of the
defined periodicity. Then, we simply calculate the scattering coefficients.

.. literalinclude:: examples/grating_tmatrixc.py
    :language: python
    :lines: 27-35

In the last step, we sum up the scattered fields at each point we want to calculate
the pressure field.

.. literalinclude:: examples/grating_tmatrixc.py
    :language: python
    :lines: 37-54

and plot the results.

.. plot:: examples/grating_tmatrixc.py


Two-dimensional arrays (in the xy plane)
========================================

As the last example, let us examine a structure that is a crystal consisting
of infinitely long cylinders in a square array in the xy plane.

.. literalinclude:: examples/crystal_tmatrixc.py
    :language: python
    :lines: 7-24

Similarly to the case of spheres and a three-dimensional lattice, we can check the
smallest singular value.

.. plot:: examples/crystal_tmatrixc.py