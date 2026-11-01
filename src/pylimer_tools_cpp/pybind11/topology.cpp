#ifndef PYBIND_ENTANGLEMENT_DETECT_H
#define PYBIND_ENTANGLEMENT_DETECT_H

#include "../topo/EntanglementDetector.h"
#include "../topo/ResistanceDistanceAnalysis.h"

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
  py::class_<AtomPairEntanglements, py::smart_holder>(
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

  py::class_<pylimer_tools::topo::ResistanceDistanceAnalysis, py::smart_holder>(
    m, "ResistanceDistanceAnalysis")
    .def(py::init<const std::vector<size_t>&,
                  const std::vector<size_t>&,
                  const std::vector<double>&,
                  size_t>(),
         py::arg("spring_from"),
         py::arg("spring_to"),
         py::arg("weights") = std::vector<double>{},
         py::arg("num_nodes") = 0)
    .def(py::init<const Eigen::SparseMatrix<double>&>(), py::arg("L"))
    .def_static(
      "create_linear_chain",
      [](size_t num_nodes) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::createLinearChain(
              num_nodes)));
      },
      py::arg("num_nodes"))
    .def_static(
      "create_ring",
      [](size_t num_nodes) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::createRing(
              num_nodes)));
      },
      py::arg("num_nodes"))
    .def_static(
      "create_complete",
      [](size_t num_nodes) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::createComplete(
              num_nodes)));
      },
      py::arg("num_nodes"))
    .def_static(
      "create_star",
      [](size_t num_arms, size_t arm_length) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::createStar(
              num_arms, arm_length)));
      },
      py::arg("num_arms"),
      py::arg("arm_length") = 1)
    .def("get_nr_of_nodes",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getNrOfNodes)
    .def("get_nr_of_edges",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getNrOfEdges)
    .def("get_total_edge_weight",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getTotalEdgeWeight)
    .def("is_connected",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::isConnected)
    .def("get_nr_of_connected_components",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           getNrOfConnectedComponents)
    .def(
      "get_connected_components",
      &pylimer_tools::topo::ResistanceDistanceAnalysis::getConnectedComponents)
    .def("get_laplacian",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getLaplacian)
    .def("compute_effective_resistance",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeEffectiveResistance,
         py::arg("i"),
         py::arg("j"))
    .def("compute_effective_resistance_batched",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeEffectiveResistanceBatched,
         py::arg("pairs"))
    .def("compute_all_pairs_effective_resistance",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeAllPairsEffectiveResistance,
         py::arg("chunk_size") = 128)
    .def("compute_moore_penrose_pseudo_inverse",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeMoorePenrosePseudoInverse,
         py::arg("chunk_size") = 128)
    .def("solve_identity_chunked",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::solveIdentityChunked,
         py::arg("chunk_size") = 128)
    .def(
      "compute_kirchhoff_index",
      &pylimer_tools::topo::ResistanceDistanceAnalysis::computeKirchhoffIndex)
    .def("compute_kirchhoff_index_spectral",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeKirchhoffIndexSpectral)
    .def("compute_total_resistance_per_node",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeTotalResistancePerNode)
    .def("compute_resistance_centrality",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeResistanceCentrality)
    .def("compute_commute_time",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::computeCommuteTime,
         py::arg("i"),
         py::arg("j"))
    .def("compute_all_pairs_commute_time",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeAllPairsCommuteTime);
};
// }

#endif
