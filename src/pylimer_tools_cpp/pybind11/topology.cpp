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
  py::class_<AtomPairEntanglements>(
    m, "AtomPairEntanglements", py::module_local(), R"pbdoc(
      A struct to store pairs of atoms that are close together and could be entanglements.
    )pbdoc")
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

    :param universe: The universe of atoms from which to sample entanglements from.
    :param nr_of_samples: The number of pairs of atoms to randomly sample.
    :param upper_cutoff: The maximum distance between atoms for a pair to be considered a potential entanglement.
    :param lower_cutoff: The minimum distance between atoms for a pair to be considered a potential entanglement.
    :param minimum_nr_of_samples: The minimum number of entanglements to be found.
    :param same_strand_cutoff: The maximum distance between atoms on the same strand for a pair to be considered a potential entanglement.
    :param seed: A seed for the random number generator.
    :param crosslinker_type: The type of crosslinker to consider when finding entanglements. Used for the splitting into strands.
    :param ignore_crosslinks: Whether to ignore crosslinks when finding entanglements. Careful: if you don't ignore them, the same-strand policy might not work correctly, since each crosslink should actually be associated with more than one strand.
    :param filter_dangling_and_soluble: Whether to filter out dangling chains and soluble crosslinks when finding entanglements.
      This means, entanglements involving an obviously (1st order) dangling or soluble chain are 
  )pbdoc",
        py::arg("universe"),
        py::arg("nr_of_samples"),
        py::arg("upper_cutoff"),
        py::arg("lower_cutoff") = 0,
        py::arg("minimum_nr_of_samples") = 0,
        py::arg("same_strand_cutoff") = 3.,
        py::arg("seed") = "",
        py::arg("crosslinker_type") = 2,
        py::arg("ignore_crosslinks") = true,
        py::arg("filter_dangling_and_soluble") = false);
};
// }

#endif
