
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_pylimer_bound_entities(py::module_ &);
void init_pylimer_bound_calc(py::module_ &);

PYBIND11_MODULE(pylimer_tools_cpp, m)
{
  init_pylimer_bound_entities(m);
  init_pylimer_bound_calc(m);
}
