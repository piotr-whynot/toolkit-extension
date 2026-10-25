# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

# Configuration file for the Sphinx documentation builder.

import tomllib  # Use 'import tomllib' if using Python 3.11+
from pathlib import Path
import os

import faulthandler
faulthandler.enable()


import os

def _exit_cleanly(app, exception):
    if exception is None:      # build succeeded: skip the crashing interpreter teardown
        os._exit(0)


# -- Project information

project = 'ACACIA_S2S_TOOLKIT'
copyright = '2025, Joshua Talib + ACACIA contributors'
author = 'Joshua Talib'

release = '1.0'
version = '1.1.4'

# -- General configuration
html_theme = "sphinx_rtd_theme"

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'myst_parser'
    ]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

nb_execution_mode = "off"

# -- Options for HTML output

html_logo = os.path.abspath("acacia_logo.png")
autodoc_mock_imports = []

# -- Options for EPUB output
epub_show_urls = 'footnote'

exclude_patterns = [
    'notebooks/*.ipynb',
]

def setup(app):
    app.connect("build-finished", _exit_cleanly)

