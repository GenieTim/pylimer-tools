
#include "version_config.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

void
init_pylimer_bound_calc(py::module_&);
void
init_pylimer_bound_entities(py::module_&);
void
init_pylimer_bound_generators(py::module_&);
void
init_pylimer_bound_readers(py::module_&);
void
init_pylimer_bound_sim(py::module_&);
void
init_pylimer_bound_topo(py::module_&);
void
init_pylimer_bound_writers(py::module_&);

PYBIND11_MODULE(pylimer_tools_cpp, m)
{
  m.doc() = R"pbdoc(
    pylimer_tools_cpp
    -----------------

    A collection of utility python functions for handling LAMMPS output and polymers in Python.

    .. autosummary::
        :toctree: _generate

    )pbdoc";

  init_pylimer_bound_readers(m);
  init_pylimer_bound_entities(m);
  init_pylimer_bound_writers(m);
  init_pylimer_bound_generators(m);
  init_pylimer_bound_topo(m);
  init_pylimer_bound_sim(m);
  init_pylimer_bound_calc(m);

  m.def(
    "version_information",
    []() {
      return "pylimer_tools, version " + std::string(__PROJECT_VERSION__) +
             "(" + std::string(__LIB_VERSION__) + "), compiled " +
             std::string(__DATE__) + " " + std::string(__TIME__) + " from " +
             std::string(__GIT_BRANCH__) + " (" + std::string(__GIT_HASH__) +
             ").";
    },
    R"pbdoc(
    Returns  a string of the the current version, incl. git hash and date of compilation.
  )pbdoc");
  m.attr("__version__") = __PROJECT_VERSION__;
}
