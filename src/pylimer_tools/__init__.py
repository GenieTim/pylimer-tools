"""
The pylimer_tools package is the non-compiled part of this library.
While the two packages (pylimer_tools and pylimer_tools_cpp) do interact and don't work standalone,
this distinction is nonetheless useful for the structure and functionality of the library,
as well as to simplify the compilation process.

In this Python part of the library,
various submodules offer functions that abstract the pylimer_tools_cpp library further,
that builds a bridge to other Python libraries, such as pandas,
or offer functionality that does not need to be as optimized for performance.
"""

import os
import sys
import pylimer_tools_cpp

sys.path.append(os.path.dirname(__file__))

__version__ = pylimer_tools_cpp.__version__
