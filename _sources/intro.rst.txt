========
Examples
========

*acoustotreams* is a program that covers various aspects of T-matrix calculations and
associated topics in acoustic scattering scenarios. *acoustotreams* adopts the functionality 
introduced, for the first time, in *treams* in the electromagnetic domain.
The functionality involves three levels: low-level functions, intermediate-level functions, 
and high-level functions and classes.

The low-level functions implement the underlying mathematical functions, that build
the foundation of T-matrix calculations. They are mainly imported from
:mod:`treams.special` and :mod:`treams.lattice`. The first one contains, e.g., 
the extensions of mathematical functions in :mod:`scipy.special`, the
second subpackage contains functions that are associated with computations in lattices.
All functions necessary in acoustic computations, e.g., the various solutions 
to the Helmholtz equation and their translation coefficients, implemented within
:mod:`acoustotreams`.

On the intermediate-level those underlying functions are combined to provide functions
as they are often needed for T-matrix calculations, e.g., the acoustic Mie coefficients 
or the expansion coefficients of scalar plane waves into spherical waves. 

The high-level functionality is more focused on the usability. We attempt to create an
useful interface to the underlying functions, that reduces redundancy and that is less
error prone than using the pure functions, while still integrating nicely with :mod:`numpy` 
functions. It consists of a combination of different classes and functions. At first,
there are the different basis sets, that can be used together with other important
parameters for the calculation, for example the embedding materials or the lattice
definitions. Then, there are the "physics-aware arrays", which keep track of these
parameters during the computation, and operators that can be applied to these arrays.
Finally, we will introduce, how these previous concepts can be applied to acoustic T-Matrices 
for scalar spherical and cylindrical solutions and to acoustic S-Matrices 
for scalar plane wave solutions.

.. toctree::
    :maxdepth: 1

    params
    physicsarray
    operators
    tmatrix
    tmatrixc
    smatrix