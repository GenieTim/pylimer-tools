Installation Guide
==================

pylimer-tools can be installed in several ways depending on your needs and system configuration.

Quick Installation (Recommended)
-------------------------------

For most users, the easiest way to install pylimer-tools is via PyPI:

.. code-block:: bash

   pip install pylimer-tools

This will install the latest stable version with pre-compiled binaries for most platforms.

Development Installation
------------------------

If you want to contribute to pylimer-tools or need the latest development features, you can install from source.

Prerequisites
~~~~~~~~~~~~~

Before installing from source, ensure you have the following development tools:

**On macOS/Linux:**

.. code-block:: bash

   # Install git and build tools
   # On macOS (using Homebrew):
   brew install git cmake bison flex python3
   
   # On Ubuntu/Debian:
   sudo apt-get update
   sudo apt-get install git build-essential cmake python3-dev bison flex

   # On CentOS/RHEL/Fedora:
   sudo yum install git gcc-c++ cmake python3-devel bison flex
   # or for newer versions:
   sudo dnf install git gcc-c++ cmake python3-devel bison flex

**On Windows:**

1. Install `Visual Studio Community <https://visualstudio.microsoft.com/vs/community/>`_ with C++ development tools
2. Install `Git for Windows <https://git-scm.com/download/win>`_
3. Install `CMake <https://cmake.org/download/>`_

Installation from Source
~~~~~~~~~~~~~~~~~~~~~~~

1. **Clone the repository:**

   .. code-block:: bash

      git clone https://github.com/GenieTim/pylimer-tools.git
      cd pylimer-tools

2. **Install in development mode:**

   .. code-block:: bash

      pip install -e .

   This installs the package in "editable" mode, so changes to the source code are immediately reflected.

3. **Verify the installation:**

   .. code-block:: python

      import pylimer_tools_cpp
      print(pylimer_tools_cpp.version_information())

Updating Your Installation
-------------------------

**For PyPI installations:**

.. code-block:: bash

   pip install --upgrade pylimer-tools

**For development installations:**

.. code-block:: bash

   cd /path/to/pylimer-tools
   git pull origin main
   pip install -e .

Getting Help
------------

If you encounter issues during installation:

1. Check the `GitHub Issues page <https://github.com/GenieTim/pylimer-tools/issues>`_
2. Review the sections above whether you have all dependencies installed
3. Create a new issue with:
   - Your operating system and Python version
   - The complete error message
   - Steps you've already tried

Virtual Environment Setup (Recommended)
--------------------------------------

For the best experience and to avoid conflicts with other packages, use a virtual environment:

.. code-block:: bash

   # Create a new virtual environment
   python -m venv pylimer_env
   
   # Activate the environment
   # On macOS/Linux:
   source pylimer_env/bin/activate
   
   # On Windows:
   pylimer_env\Scripts\activate
   
   # Install pylimer-tools
   pip install pylimer-tools
   
   # When done, deactivate
   deactivate

With conda/mamba:

.. code-block:: bash

   # Create a new conda environment
   conda create -n pylimer_env python=3.13
   conda activate pylimer_env
   
   # Install pylimer-tools
   pip install pylimer-tools
