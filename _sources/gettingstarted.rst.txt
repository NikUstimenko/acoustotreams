.. highlight:: console

===============
Getting started
===============

Installation
============

To install the package with pip, use ::

   pip install acoutreams


How to use acoutreams
=====================

Import *acoutreams*, create acoustic T-matrices and start calculating.

.. doctest::

   >>> import acoutreams
   >>> tm = acoutreams.AcousticTMatrix.sphere(3, 600, 0.01, [(1.29, 343), (1500, 1000)])
   >>> f"{tm.xs_ext_avg:.6f}"
   '0.000941'

More detailed examples are given in :doc:`intro`.