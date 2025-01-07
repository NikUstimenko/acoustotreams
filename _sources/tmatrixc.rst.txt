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

Cylindrical T-matrices for one-dimensional arrays of spherical T-matrices
=========================================================================

One-dimensional arrays (along the x axis)
=========================================

Two-dimensional arrays (in the xy plane)
=========================================