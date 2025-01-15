.. acoustotreams documentation master file, created by
   sphinx-quickstart on Thu Dec 19 15:58:08 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to acoustotreams's documentation!
======================================

.. toctree::
   :maxdepth: 1

   gettingstarted
   intro
   theory
   dev
   about
   

The package **acoustotreams** adopts the framework of the package **treams** for
electromagnetic scattering computations to simplify computations of the
scattering of acoustic waves in acoustic metamaterials. The supported geometries 
include single particles as well as finite and periodically infinite arrangements. 
All methods are suitable for the use of lossy materials. The periodic systems can have 
one-, two-, or three-dimensional lattices. The lattice computations are performed 
by the functions imported from :mod:`treams.lattice` which accelerates them 
by converting the occurring slowly converging summations to exponentially fast convergent series. 
From :mod:`treams.special` the mathematical functions, which are typically necessary 
in T-Matrix method computations, are also imported.

To accommodate the periodic structures of different dimensionalities, three types of
solutions to the scalar Helmholtz equation are employed: plane waves, cylindrical
waves, and spherical waves. For each of those solution sets, the typical manipulations,
e.g. translations and rotations, are implemented, as well as transformations between
them.

Finally, three classes are the main point of interaction for the user. They allow access
to the underlying functions operating directly on the spherical and cylindrical
acoustic T-matrices or the acoustic S-matrices based on the plane wave solutions.

.. todo:: clean up intro

Features
========

* T-matrix calculations using a spherical or cylindrical wave basis set
* Scattering from clusters of particles
* Scattering from particles and clusters arranged in 3d-, 2d-, and 1d-lattices
* Calculation of sound propagation in stratified media
* Band calculation in crystal structures

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
