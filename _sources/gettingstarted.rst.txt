.. highlight:: console

===============
Getting started
===============

Installation
============

To install the package with pip, use ::

   pip install acoustotreams


How to use acoustotreams
========================

Import *acoustotreams*, create acoustic T-matrices and start calculating.

.. doctest::

   >>> import acoustotreams
   >>> tm = acoustotreams.AcousticTMatrix.sphere(2, 600, 0.01, [(1.3, 343), (998, 1497)])
   >>> f"{tm.xs_ext_avg:.6f}"
   '0.001011'

More detailed examples are given in :doc:`intro`.