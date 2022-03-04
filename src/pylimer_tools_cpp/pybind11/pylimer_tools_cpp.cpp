
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_pylimer_bound_entities(py::module_ &);
void init_pylimer_bound_calc(py::module_ &);
void init_pylimer_bound_readers(py::module_ &);
void init_pylimer_bound_writers(py::module_ &);
void init_pylimer_bound_generators(py::module_ &);

PYBIND11_MODULE(pylimer_tools_cpp, m) {
  m.doc() = R"pbdoc(
    PylimerTools Cpp
    -----------------

    A collection of utility python functions for handling LAMMPS output and polymers in Python.

    .. autosummary::
        :toctree: _generate

    )pbdoc";

  init_pylimer_bound_entities(m);
  init_pylimer_bound_readers(m);
  init_pylimer_bound_writers(m);
  init_pylimer_bound_generators(m);
  init_pylimer_bound_calc(m);
}
