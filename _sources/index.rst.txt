.. acoustotreams documentation master file, created by
   sphinx-quickstart on Thu Dec 19 15:58:08 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to acoustotreams's documentation!
=========================================

.. toctree::
   :maxdepth: 1

   gettingstarted
   intro
   theory
   acoustotreams
   dev
   about
   

The package **acoustotreams** adopts the framework of the package **treams** to 
simplify computations of the scattering of acoustic waves in acoustic metamaterials. 
The supported geometries include single particles as well as finite and periodically infinite arrangements. 
All methods are suitable for the use of lossy materials. The periodic systems can have 
one-, two-, or three-dimensional lattices. The lattice computations are performed 
by the functions imported from :mod:`treams.lattice`, which accelerates them 
by converting the occurring slowly converging sums to exponentially fast convergent series
using the Ewald summation technique. From :mod:`acoustotreams.special` the mathematical functions, 
which are typically necessary in T-matrix computations, are also imported.

To accommodate acoustic scatterers of different symmetry and also the periodic structures 
of different dimensionalities, three types of solutions to the scalar Helmholtz equation 
are employed: plane waves, cylindrical waves, and spherical waves. For each of those solution sets, 
the typical manipulations, e.g., translations and rotations, are implemented, as well as 
transformations between them.

Finally, three classes are the main point of interaction for the user. They allow for access
to the underlying functions that directly transorm acoustic T-matrices in the spherical- or 
cylindrical-wave basis and acoustic S-matrices in the plane-wave basis.

.. todo:: clean up intro

Features
========

* T-matrix-based computations using a spherical- or cylindrical-wave basis set
* Scattering by clusters of scatterers
* Scattering by individiual scatterers or clusters thereof arranged in 3d-, 2d-, and 1d-lattices
* Simulation of propagation of acoustic fields in stratified media
* Band calculation for crystal structures

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
