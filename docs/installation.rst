Installation
============

You can use PyPip to install this package:

.. code:: bash

  python -m pip install pylimer_tools

Note that, depending on your OS and Python environment, you may require an installation of the igraph library. 

Error: igraph not found
-----------------------

An OS-independent installation of the igraph library can be done using `vcpkg`_ like so:

.. code:: bash

  vcpkg install igraph

For CMake to find the installed library, be sure to expose the root folder of `vcpkg` as the environment variable `VCPKG_ROOT`.

On MacOS, you can install the igraph library using brew:

.. code:: bash

  brew install igraph

.. _vcpkg: https://github.com/microsoft/vcpkg
