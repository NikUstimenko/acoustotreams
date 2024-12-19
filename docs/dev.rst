. highlight:: console

============================
Development and contributing
============================

Setting up the development environment
======================================

1.
Clone the repository with ::

   git clone git@github.com:NikUstimenko/acoutreams.git

or ::

   git clone https://github.com/NikUstimenko/acoutreams.git

and enter the directory ::

   cd acoutreams

2.
This step is optional, but you might want to create an environment where acoutreams and the
development tools are created ::

   conda env create --name acoutreams-dev python

Activate the environment with::

   conda activate acoutreams-dev

3.
Setup the package with ::

   pip install -e .

This last step makes the program available in the environment independently of the
current folder. This is especially necessary for correctly building the documentation.


Building the documentation
==========================


To install the required packages for testing use ::

   pip install acoutreams[docs]

The documentation is built with sphinx by using ::

   sphinx-build -b html docs docs/_build/html

from the root directory of the package to build the documentation as html pages.
Some figures are automatically created, which requires matplotlib.

The doctests can be run with ::

   sphinx-build -b doctest docs docs/_build/doctest
