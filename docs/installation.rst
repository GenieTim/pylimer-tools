Installation
============

You can use PyPip to install this package:

.. code:: bash

  python -m pip install pylimer_tools

**Note**: depending on your OS and Python environment, you may require an installation of some dependencies.
In this case, a custom compilation could be a workaroud.
Checkout the following sections for your system. 

MacOS / Unix systems
---------------------

A complete script to install this library including all dependencies could look like this:

.. code:: bash

  # install vcpkg
  git clone https://github.com/microsoft/vcpkg
  ./vcpkg/bootstrap-vcpkg.sh

  # install dependencies
  ./vcpkg/vcpkg install igraph
  ./vcpkg/vcpkg install boost

  # indicate that vcpkg should be used as toolchain
  export VCPKG_ROOT="$(pwd)/vcpkg"

  # now, the compilation should work
  git clone https://github.com/GenieTim/pylimer-tools pylimer_tools
  python -m pip install ./pylimer_tools

To update the custom installation, follow the steps below:

.. code:: bash
  
  cd pylimer_tools
  git pull
  python -m pip install .

Windows
--------

Please make sure to have a current version of `Visual Studio`_ as well as `git`_ installed before following the next steps.
A complete PowerShell script to install this library including all dependencies could look like this:

.. code:: bash

  # install vcpkg
  git clone https://github.com/microsoft/vcpkg
  .\vcpkg\bootstrap-vcpkg.sh

  # install dependencies
  .\vcpkg\vcpkg install igraph
  .\vcpkg\vcpkg install boost

  # indicate that vcpkg should be used for linking
  $env:VCPKG_ROOT="$(pwd)/vcpkg"

  # now, the compilation should work
  git clone https://github.com/GenieTim/pylimer-tools pylimer_tools
  python -m pip install .\pylimer_tools

To update the custom installation, follow the steps below:

.. code:: bash

  $env:VCPKG_ROOT="$(pwd)/vcpkg"
  cd pylimer_tools
  git pull
  python -m pip install .

Frequent Errors
---------------

igraph not found
~~~~~~~~~~~~~~~~~

This means, the library code:`igraph` that is used by this library was not found. 

An OS-independent installation of the igraph library can be done using `vcpkg`_ like so:

.. code:: bash

  vcpkg install igraph
  vcpkg install boost

For CMake to find the installed library, be sure to expose the root folder of `vcpkg` as the environment variable `VCPKG_ROOT`.
Note that on Windows, vcpkg has multiple triplets that it chooses from when finding the packages. 
It might be possible that you have igraph/boost installed for one triplet, but another one is searched for.
Specify `./vcpkg.exe install igraph --triplet=x64-windows` or `./vcpkg.exe install igraph --triplet=x86-windows` to install the other one too.

On MacOS, you can install the igraph library using brew:

.. code:: bash

  brew install igraph




.. _vcpkg: https://github.com/microsoft/vcpkg
.. _git: https://www.git-scm.com/
.. _Visual Studio: https://visualstudio.microsoft.com/
