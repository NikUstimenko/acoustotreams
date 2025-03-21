.. highlight:: python

.. only:: builder_html

===================
Acoustic T-Matrices
===================

.. contents:: Table of contents
   :local:

One of the main objects for acoustic scattering calculations within *acoustotreams* are
acoustic T-matrices. They describe the scattering response of an object by encoding the linear
relationship between incident and scattered fields. These fields are expanded using
the scalar spherical waves.

The acoustic T-matrices of spheres can be obtained analytically. For more complicated shapes
numerical methods are necessary to compute the acoustic T-matrix. Once the T-matrix of a single
object is known, the acoustic interaction between particle within clusters can be
calculated efficiently. Such clusters can be analyzed in their local description, where
the field expansions are centered at each particle of the cluster, or in a global
description treating the whole cluster as a single object.

*acoustotreams* is particularly aimed at analyzing scattering within lattices. These lattices
can be periodic in one, two, or all three spatial dimensions. The unit cell of those
lattices can consist of an arbitrary number of objects described by the acoustic T-matrix.

Spheres
=======

It's possible to calculate the acoustic T-matrix of a single sphere with the
method :meth:`~acoustotreams.AcousticTMatrix.sphere`. We start by defining the relevant parameters for
our calculation and creating the acoustic T-matrices themselves.

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 7-12

Here we define a material using the class :meth:`~acoustotreams.AcousticMaterial`.
To create an instance of this class, we pass three arguments: mass density,
longitudinal speed of sound, and transverse speed of sound. The acoustic
material parameters can be complex. The imaginary part of the density must be
positive, whereas the imaginary part of the speeds negative. 

Now, we can easily access quantities like the scattering and extinction cross
sections

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 14-15

From the parameter ``lmax = 10`` we see that the acoustic T-matrix is calculated up to the tenth
multipolar order. To restrict the T-matrix in the monopole-dipole approximation, we can
select a basis containing only those multipoles.

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 17-20

Now, we can look at the results by plotting them and observe, unsurprisingly, that for
larger frequencies the monopole-dipole approximation is not giving an accurate result.
Further, we visualize the total acoustic field intensity at frequency 163 kHz.

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 22-40

We select the T-matrix and illuminate it with a plane wave. Next, we set up 
the x and z coordinates and define the auxiliary function to compute
the intensity at a single valid point. We can calculate the intensity
of the fields as a superposition of incident and scattered fields.
Here, we also accelerate the computations using paralyzation by
:mod:`Parallel` in :mod:`joblib`. 

Finally, we compute the radiation pattern of the sphere at the same frequency
as a function of the polar angle :math:`\theta`. The red arrow indicates 
the direction of incidence of the plane wave.

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 42-47

.. plot:: examples/sphere.py


Clusters
========

Multi-scattering calculations in a cluster of particles is a typical application of the
T-matrix method. We first construct an object from two spheres. Using the definition 
of the relevant parameters

.. literalinclude:: examples/cluster.py
   :language: python
   :lines: 7-17

we can simply first create the spheres and put them together in a cluster, where we immediately calculate 
the interaction.

.. literalinclude:: examples/cluster.py
   :language: python
   :lines: 19-20

Then, we can illuminate with a plane wave and get the scattered field coefficients and
the scattering and extinction cross sections for that particular illumination.

.. literalinclude:: examples/cluster.py
   :language: python
   :lines: 22-24

Finally, with few lines similar to the plotting of the field intensity of a single
sphere we can obtain the fields outside of the sphere.

.. literalinclude:: examples/cluster.py
   :language: python
   :lines: 26-43

Up to here, we did all calculations for the cluster in the local basis. By expanding
the incident and scattered fields in a basis with a single origin we can describe the
same object. Often, a larger number of multipoles is needed to do so and some
information on fields between the particles is lost. However, the description in a global
basis can be more efficient in terms of matrix size.

.. literalinclude:: examples/cluster.py
   :language: python
   :lines: 76-77

A comparison of the calculated near-fields and the cross sections show good agreement
between the results of both, local and global, T-matrices.

In the last figure, the T-matrix is rotated by 90 degrees about the y axis and the
illumination is set accordingly to be a plane wave propagating in the x direction, 
such that the whole system remains the same. It shows how the rotate operator 
produces consistent results.

.. literalinclude:: examples/cluster.py
   :language: python
   :lines: 117-118

.. plot:: examples/cluster.py

Clusters (Born approximations)
==============================

One-dimensional arrays (along z)
================================

Next, we turn to systems that are periodic in the z direction. We calculate the
scattering from an array of spheres. Intentionally, we choose a unit cell with two
different spheres that overlap along the z direction, but are not placed exactly along 
the same line. This is the most general case for the implemented lattice sums. 
After the common setup of the parameters, we simply create a cluster in a local basis.

.. literalinclude:: examples/chain.py
   :language: python
   :lines: 7-15

This time we let them interact specifying a one-dimensional lattice, so that the spheres
form a chain.

.. literalinclude:: examples/chain.py
   :language: python
   :lines: 17-18

Next, we choose set the illumination to be propagating along the x axis.
The z component of the wave vector of the plane wave has to match to the wave vector
component of the lattice interaction, i.e., the Bloch wave vector.

.. literalinclude:: examples/chain.py
   :language: python
   :lines: 20-21

There is an efficient way to calculate the acoustic response, especially in the far-field,
using cylindrical acoustic T-matrices. That will be introduced in :doc:`tmatrixc`. Here, we will
stay in the expression of the fields as scalar spherical waves. This allows the
calculation of the fields in the domain between the spheres. To compute them accurately, we
expand the scattered field coefficients in the whole lattice in monopole approximation at each point
we want to calculate the pressure field.

.. literalinclude:: examples/chain.py
   :language: python
   :lines: 23-40

To calculate the velocity field, we again expand the scattered field coefficients, 
however, in the dipole approximation.

.. literalinclude:: examples/chain.py
   :language: python
   :lines: 65-80

.. plot:: examples/chain.py

Two-dimensional arrays (in the xy plane)
========================================

The case of periodicity in two directions is similar to the case of the previous section
with one-dimensional periodicity. Here, by convention the array has to be in the
xy plane.

.. literalinclude:: examples/grid.py
   :language: python
   :lines: 7-42

With a few changes we computed the fields in a square array of the same spheres as in the
previous examples. Most importantly we changed the value of the variable :code:`lattice`
to an instance of a two-dimensional :class:`~acoustotreams.Lattice` and set :code:`kpar`
accordingly. Most other changes are just resulting from the change of the coordinate
system.

.. plot:: examples/grid.py

Three-dimensional arrays
========================

In a three-dimensional lattice we are mostly concerned with finding eigenmodes of a
crystal. We want to restrict the example to calculating modes at the gamma point in
reciprocal space. The calculated system consists of a single sphere in a cubic lattice.
In our very crude analysis, we blindly select the lowest singular value of the lattice
interaction matrix. Besides the mode when the frequency tends to zero, there are two
additional modes at higher frequencies in the chosen range.

.. literalinclude:: examples/crystal.py
   :language: python
   :lines: 7-18

.. plot:: examples/crystal.py