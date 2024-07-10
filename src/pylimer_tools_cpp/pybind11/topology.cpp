#ifndef PYBIND_ENTANGLEMENT_DETECT_H
#define PYBIND_ENTANGLEMENT_DETECT_H

#include "../topo/EntanglementDetector.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;

using namespace pylimer_tools::topo::entanglement_detection;
// namespace pylimer_tools::calc {
void
init_pylimer_bound_topo(py::module_& m)
{
  py::class_<AtomPairEntanglements>(m, "AtomPairEntanglements")
    .def(py::init<>(), "Get an instance of this struct")
    .def_readwrite(
      "pairs_of_atoms", &AtomPairEntanglements::pairsOfAtoms, R"pbdoc(
      A list of pairs of atom ids that are close together and could be entanglements
    )pbdoc")
    .def_readwrite("pair_of_atom",
                   &AtomPairEntanglements::pairOfAtom,
                   R"pbdoc(
      An index in the pairs_of_atoms if the atom is part of a pair, -1 else.
    )pbdoc");

  m.def("randomly_sample_entanglements",
        &randomlyFindEntanglements,
        R"pbdoc(
    Randomly find pairs of atoms that are close together and could be
    entanglements
  )pbdoc",
        py::arg("universe"),
        py::arg("nr_of_samples"),
        py::arg("cutoff"),
        py::arg("minimum_nr_of_samples") = 0,
        py::arg("same_strand_cutoff") = 3.,
        py::arg("seed") = "",
        py::arg("crosslinker_type") = 2,
        py::arg("ignore_crosslinks") = true);
}
// }

#endif
