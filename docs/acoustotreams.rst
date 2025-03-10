=========
Reference
=========

.. automodule:: acoustotreams
   :no-members:
   :no-inherited-members:
   :no-special-members:


Modules
=======

These modules provide basic functionality for transformations within one basis set, i.e.
one module, like translations and rotations as well as transformations among them.
The functions in there provide an intermediate stage between the purely mathematical
functions, either implemented within :mod:`acoustotreams` or imported from 
:mod:`treams.special` and :mod:`treams.lattice`, and the
higher-level classes and functions.

.. autosummary::
   :toctree: generated/
   :template: bare-module

   ~acoustotreams.ssw
   ~acoustotreams.scw
   ~acoustotreams.spw

Finally, a module for calculating scattering coefficients.

.. autosummary::
   :toctree: generated/

   ~acoustotreams.coeffs
