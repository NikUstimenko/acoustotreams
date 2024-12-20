.. highlight:: python

.. only:: builder_html

==========
Acoustic T-Matrices
==========

.. contents:: Table of contents
   :local:

One of the main objects for acoustic scattering calculations within *acoutreams* are
acoustic T-matrices. They describe the scattering response of an object by encoding the linear
relationship between incident and scattered fields. These fields are expanded using
the scalar spherical waves.

The acoustic T-matrices of spheres can be obtained analytically. For more complicated shapes
numerical methods are necessary to compute the acoustic T-matrix. Once the T-matrix of a single
object is known, the acoustic interaction between particle within clusters can be
calculated efficiently. Such clusters can be analyzed in their local description, where
the field expansions are centered at each particle of the cluster, or in a global
description treating the whole cluster as a single object.

*acoutreams* is particularly aimed at analyzing scattering within lattices. These lattices
can be periodic in one, two, or all three spatial dimensions. The unit cell of those
lattices can consist of an arbitrary number of objects described by the acoustic T-matrix.

Spheres
=======

It's possible to calculate the acoustic T-matrix of a single sphere with the
method :meth:`~acoutreams.AcousticTMatrix.sphere`. We start by defining the relevant parameters for
our calculation and creating the acoustic T-matrices themselves.

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 7-11

Now, we can easily access quantities like the scattering and extinction cross
sections

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 13-14

From the parameter ``lmax = 10`` we see that the acoustic T-matrix is calculated up to the tenth
multipolar order. To restrict the T-matrix in the monopole-dipole approximation, we can
select a basis containing only those multipoles.

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 16-19

Now, we can look at the results by plotting them and observe, unsurprisingly, that for
larger frequencies the monopole-dipole approximation is not giving an accurate result.
Further, we visualize the total acoustic field intensity at frequency 163 kHz.

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 21-39

We select the T-matrix and illuminate it with a plane wave. Next, we set up 
the x and z coordinates and define the auxiliary function to compute
the intensity at a single valid point. We can calculate the intensity
of the fields as a superposition of incident and scattered fields.
Here we also accelerate the computations using paralyzation in :mod:`joblib`. 

Finally, we compute the radiation pattern of the sphere at the same frequency.

.. literalinclude:: examples/sphere.py
   :language: python
   :lines: 41-49

.. plot:: examples/sphere.py