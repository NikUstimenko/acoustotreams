"""Sphinx config."""
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from importlib.metadata import version

import numpy as np

# sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = "acoustotreams"
copyright = "2024, Nikita Ustimenko"
author = "Nikita Ustimenko"
release = version("acoustotreams")
version = ".".join(release.split(".")[:3])


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "matplotlib.sphinxext.plot_directive",
]

autosectionlabel_prefix_document = True
doctest_global_setup = """
import numpy as np
from acoustotreams import *

np.set_printoptions(precision=3, suppress=True)
"""
intersphinx_mapping = {
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}
todo_include_todos = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pyramid"
html_css_files = ["custom.css"]
html_sidebars = {
    "**": ["globaltoc.html", "relations.html", "sourcelink.html", "searchbox.html"]
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]