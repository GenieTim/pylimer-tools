import os
import sys

# make sure we use the "live" pylimer tools, not the installed one
# but only the Python one, for the cpp one, we need the copiled/installed
import pylimer_tools_cpp

pylimer_tools_cpp_path = pylimer_tools_cpp.__file__
sys.path.insert(0, os.path.dirname(__file__) + "/../../src")
sys.path.insert(0, os.path.dirname(pylimer_tools_cpp_path))
# sys.path.insert(0, os.path.dirname(__file__) + "/../../src/pylimer_tools/calc")
# sys.path.insert(0, os.path.dirname(__file__) + "/../../src/pylimer_tools/io")
# sys.path.insert(0, os.path.dirname(__file__) +
#                 "/../../src/pylimer_tools/utils")

sys.path.insert(0, os.path.dirname(__file__))
