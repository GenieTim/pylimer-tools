# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys
import sysconfig
import warnings
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from datetime import datetime

from intersphinx_registry import get_intersphinx_mapping

import pylimer_tools_cpp

# -- Project information -----------------------------------------------------

project = "pylimer-tools"
copyright = "2021-" + datetime.now().strftime("%Y") + ", Tim Bernhard"
author = "Tim Bernhard"

# The full version, including alpha/beta/rc tags
version = pylimer_tools_cpp.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_automodapi.automodapi",
    "sphinx_copybutton",
    "sphinx_gallery.gen_gallery",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.githubpages",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinxext.opengraph",
]

# Configuration for specific extensions
# napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# opengraph settings
ogp_site_url = "https://genietim.github.io/pylimer-tools"
ogp_custom_meta_tags = [
    '<link rel="icon" type="image/png" href="/_static/favicon/favicon-96x96.png" sizes="96x96" />',
    '<link rel="icon" type="image/svg+xml" href="/_static/favicon/favicon.svg" />',
    '<link rel="shortcut icon" href="/_static/favicon/favicon.ico" />',
    '<link rel="apple-touch-icon" sizes="180x180" href="/_static/favicon/apple-touch-icon.png" />',
    '<meta name="apple-mobile-web-app-title" content="pylimer-tools" />',
    '<link rel="manifest" href="/_static/favicon/site.webmanifest" />',
]

# sphinx-gallery settings
sphinx_gallery_conf = {
    "examples_dirs": "../examples",  # path to your example scripts
    "gallery_dirs": "auto_examples",  # path to where to save gallery generated output
    "show_signature": True,
}

intersphinx_mapping = get_intersphinx_mapping(
    packages={
        "matplotlib",
        "numpy",
        "pandas",
        "python",
    },
)

# To avoid showing methods and attributes of classes multiple times.
numpydoc_show_class_members = False

# Enable autosummary
autosummary_generate = True
autoclass_content = "both"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = [".rst", ".md"]

# The master toctree document.
master_doc = "index"

# Mathjax options
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"
# mathjax3_config = {
#     "tex": {
#         # "inlineMath": [['$', '$'], ['\\(', '\\)']]
#     },
#     "extensions": ["jsMath2jax.js"],
#     "jax": ["input/TeX"],
# }
# mathjax_options = {"async": "async"}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

# Theme-specific options
html_theme_options = {
    "sidebar_hide_name": True,
    "navigation_with_keys": True,
    "top_of_page_buttons": ["view", "edit"],
    "source_repository": "https://github.com/GenieTim/pylimer-tools",
    "source_branch": "main",
    "source_directory": "docs/",
    "logo": {
        "text": "pylimer-tools",
    },
    "secondary_sidebar_items": ["page-toc", "sg_download_links", "sg_launcher_links"],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Custom CSS files
html_css_files = [
    "css/custom.css",
]

# The name for this set of Sphinx documents.
html_title = f"{project} v{version}"

# A shorter title for the navigation bar.
html_short_title = f"{project}"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = None

# If false, no module index is generated.
html_domain_indices = True

# If false, no index is generated.
html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',
    # Latex figure (float) alignment
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
# latex_documents = [
#     (
#         master_doc,
#         "python_example.tex",
#         "python_example Documentation",
#         "Sylvain Corlay",
#         "manual",
#     ),
# ]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
# latex_use_parts = False

# If true, show page references after internal links.
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
# latex_domain_indices = True
