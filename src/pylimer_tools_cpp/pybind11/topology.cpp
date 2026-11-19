#ifndef PYBIND_ENTANGLEMENT_DETECT_H
#define PYBIND_ENTANGLEMENT_DETECT_H

#include "../entities/Universe.h"
#include "../topo/EntanglementDetector.h"
#include "../topo/ResistanceDistanceAnalysis.h"
#include "../utils/CrosslinkerUniverseStats.h"

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;

using namespace pylimer_tools::topo::entanglement_detection;

// namespace pylimer_tools::calc {
void
init_pylimer_bound_topo(py::module_& m)
{
  py::native_enum<EntanglementFilteringMode>(
    m,
    "EntanglementFilteringMode",
    "enum.IntEnum",
    "Enum representing the filtering mode for entanglement sampling.")
    .value("NONE",
           EntanglementFilteringMode::NONE,
           "No filtering: all entanglements retained")
    .value("OMIT_INACTIVE",
           EntanglementFilteringMode::OMIT_INACTIVE,
           "Melt density mode: sample from all strands, omit from output if "
           "touching dangling/soluble strands")
    .value("PRUNE_UPFRONT",
           EntanglementFilteringMode::PRUNE_UPFRONT,
           "Active-only mode: topologically peel leaves upfront, sample "
           "strictly between active strands")
    .finalize();

  m.def("peel_topological_leaves",
        &peelTopologicalLeaves,
        R"pbdoc(
    Peel 1-core dangling trees and subgraphs iteratively from the CrosslinkerUniverse topology.
    Returns a boolean vector indicating for each strand whether it is ignored (pruned).

    :param crosslinker_universe: The CrosslinkerUniverse representation.
    :return: Boolean vector of ignored/pruned strands.
  )pbdoc",
        py::arg("crosslinker_universe"));

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

  py::class_<StrandPairEntanglements, py::smart_holder>(
    m, "StrandPairEntanglements", py::module_local(), R"pbdoc(
      A struct to store sampled entanglements on strands (without requiring an explicit Universe of all beads).
    )pbdoc")
    .def(py::init<>(), "Get an instance of this struct")
    .def_readwrite(
      "strand_idx_1", &StrandPairEntanglements::strandIdx1, R"pbdoc(
      First strand index of each entanglement pair.
    )pbdoc")
    .def_readwrite(
      "strand_idx_2", &StrandPairEntanglements::strandIdx2, R"pbdoc(
      Second strand index of each entanglement pair.
    )pbdoc")
    .def_readwrite("alpha_1", &StrandPairEntanglements::alpha1, R"pbdoc(
      Contour fraction along the first strand in (0, 1).
    )pbdoc")
    .def_readwrite("alpha_2", &StrandPairEntanglements::alpha2, R"pbdoc(
      Contour fraction along the second strand in (0, 1).
    )pbdoc")
    .def_readwrite("bead_idx_1", &StrandPairEntanglements::beadIdx1, R"pbdoc(
      Bead index along the first strand in [0, beads_in_strand - 1].
    )pbdoc")
    .def_readwrite("bead_idx_2", &StrandPairEntanglements::beadIdx2, R"pbdoc(
      Bead index along the second strand in [0, beads_in_strand - 1].
    )pbdoc")
    .def_readwrite("x", &StrandPairEntanglements::x, R"pbdoc(
      X-coordinates of the sampled entanglement midpoints.
    )pbdoc")
    .def_readwrite("y", &StrandPairEntanglements::y, R"pbdoc(
      Y-coordinates of the sampled entanglement midpoints.
    )pbdoc")
    .def_readwrite("z", &StrandPairEntanglements::z, R"pbdoc(
      Z-coordinates of the sampled entanglement midpoints.
    )pbdoc")
    .def_readwrite("position_1", &StrandPairEntanglements::position1, R"pbdoc(
      Sampled 3D positions of the first entangled beads.
    )pbdoc")
    .def_readwrite("position_2", &StrandPairEntanglements::position2, R"pbdoc(
      Sampled 3D positions of the second entangled beads.
    )pbdoc")
    .def("size", &StrandPairEntanglements::size)
    .def("__len__", &StrandPairEntanglements::size)
    .def("to_atom_pair_entanglements",
         &StrandPairEntanglements::toAtomPairEntanglements,
         py::arg("crosslinker_universe"),
         R"pbdoc(
      Convert the strand pair entanglements to AtomPairEntanglements using the deterministic atom indexing of CrosslinkerUniverse.
    )pbdoc");

  m.def("sample_strand_beads_on_demand",
        &sampleStrandBeadsOnDemand,
        R"pbdoc(
    Deterministically realize the internal bead coordinates of a single strand
    in a CrosslinkerUniverse using a per-strand PRNG seed.

    :param box: The simulation box.
    :param crosslinker_universe: The CrosslinkerUniverse of crosslinkers and strands.
    :param strand_idx: Index of the strand in `crosslinker_universe`.
    :param strand_seed: Deterministic per-strand PRNG seed.
    :param n_mc_steps: Number of Monte Carlo equilibration steps (default: 0).
    :return: Vector of length `3 * beads_in_strand[strand_idx]` with bead coordinates.
  )pbdoc",
        py::arg("box"),
        py::arg("crosslinker_universe"),
        py::arg("strand_idx"),
        py::arg("strand_seed"),
        py::arg("n_mc_steps") = 0);

  m.def("randomly_sample_entanglements_on_demand",
        py::overload_cast<const pe::Box&,
                          const pylimer_tools::sim::mehp::CrosslinkerUniverse&,
                          const size_t,
                          const double,
                          const double,
                          const size_t,
                          const double,
                          const std::string&,
                          EntanglementFilteringMode,
                          int>(&randomlyFindEntanglementsOnDemand),
        R"pbdoc(
    Sample pairs of strand beads that could form entanglements on demand from a CrosslinkerUniverse
    using per-strand bounding boxes and on-demand Brownian bridge bead resampling.
    Only strand seeds and bounding boxes are stored, achieving O(N_strands) memory scaling
    while sampling from the exact discrete multi-bead chain conformations.

    :param box: The simulation box for periodic boundary conditions.
    :param crosslinker_universe: The CrosslinkerUniverse of crosslinkers and strands.
    :param nr_of_samples: The number of entanglement pairs to randomly sample.
    :param upper_cutoff: The maximum distance between two beads to be considered an entanglement.
    :param lower_cutoff: The minimum distance between two beads to be considered an entanglement.
    :param minimum_nr_of_samples: The minimum number of entanglements to be found.
    :param same_strand_cutoff: The minimum contour separation (in beads) between two beads on the same strand.
    :param seed: A seed for the random number generator.
    :param filtering_mode: Filtering mode for entanglements (default: EntanglementFilteringMode.OMIT_INACTIVE).
    :param n_mc_steps: Number of MC equilibration steps per chain (default: 0).
  )pbdoc",
        py::arg("box"),
        py::arg("crosslinker_universe"),
        py::arg("nr_of_samples"),
        py::arg("upper_cutoff"),
        py::arg("lower_cutoff") = 0.,
        py::arg("minimum_nr_of_samples") = 0,
        py::arg("same_strand_cutoff") = 3.,
        py::arg("seed") = "",
        py::arg_v("filtering_mode",
                  EntanglementFilteringMode::OMIT_INACTIVE,
                  "EntanglementFilteringMode.OMIT_INACTIVE"),
        py::arg("n_mc_steps") = 0);

  m.def(
    "randomly_sample_entanglements_on_demand",
    [](const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse,
       const pe::Box& box,
       const size_t nrOfSamples,
       const double upperCutoff,
       const double lowerCutoff,
       const size_t minimumNrOfSamples,
       const double sameStrandCutoff,
       const std::string& seed,
       const EntanglementFilteringMode filteringMode,
       const int nMcSteps) {
      return randomlyFindEntanglementsOnDemand(box,
                                               crosslinkerUniverse,
                                               nrOfSamples,
                                               upperCutoff,
                                               lowerCutoff,
                                               minimumNrOfSamples,
                                               sameStrandCutoff,
                                               seed,
                                               filteringMode,
                                               nMcSteps);
    },
    py::arg("crosslinker_universe"),
    py::arg("box"),
    py::arg("nr_of_samples"),
    py::arg("upper_cutoff"),
    py::arg("lower_cutoff") = 0.,
    py::arg("minimum_nr_of_samples") = 0,
    py::arg("same_strand_cutoff") = 3.,
    py::arg("seed") = "",
    py::arg_v("filtering_mode",
              EntanglementFilteringMode::OMIT_INACTIVE,
              "EntanglementFilteringMode.OMIT_INACTIVE"),
    py::arg("n_mc_steps") = 0);

  m.def(
    "randomly_sample_entanglements",
    py::overload_cast<const pylimer_tools::entities::Universe&,
                      const size_t,
                      const double,
                      const double,
                      const size_t,
                      const double,
                      const std::string&,
                      const int,
                      const bool,
                      EntanglementFilteringMode>(&randomlyFindEntanglements),
    R"pbdoc(
    Randomly find pairs of atoms that are close together and could be entanglements.

    :param universe: The universe of atoms from which to sample entanglements from.
    :param nr_of_samples: The number of pairs of atoms to randomly sample.
    :param upper_cutoff: The maximum distance between atoms for a pair to be considered a potential entanglement.
    :param lower_cutoff: The minimum distance between atoms for a pair to be considered a potential entanglement.
    :param minimum_nr_of_samples: The minimum number of entanglements to be found.
    :param same_strand_cutoff: The minimum shortest path distance on the graph between atoms for an entanglement.
    :param seed: A seed for the random number generator.
    :param crosslinker_type: Atom type ID representing crosslinkers.
    :param ignore_crosslinks: Whether to exclude crosslink atoms from forming entanglements.
    :param filtering_mode: Filtering mode for entanglements (default: EntanglementFilteringMode.NONE).
  )pbdoc",
    py::arg("universe"),
    py::arg("nr_of_samples"),
    py::arg("upper_cutoff"),
    py::arg("lower_cutoff") = 0.0,
    py::arg("minimum_nr_of_samples") = 0,
    py::arg("same_strand_cutoff") = 3.0,
    py::arg("seed") = "",
    py::arg("crosslinker_type") = 2,
    py::arg("ignore_crosslinks") = true,
    py::arg_v("filtering_mode",
              EntanglementFilteringMode::NONE,
              "EntanglementFilteringMode.NONE"));

  m.def("randomly_sample_entanglements",
        py::overload_cast<const pe::Box&,
                          const pylimer_tools::sim::mehp::CrosslinkerUniverse&,
                          const size_t,
                          const double,
                          const double,
                          const size_t,
                          const double,
                          const std::string&,
                          EntanglementFilteringMode,
                          int>(&randomlyFindEntanglementsOnDemand),
        py::arg("box"),
        py::arg("crosslinker_universe"),
        py::arg("nr_of_samples"),
        py::arg("upper_cutoff"),
        py::arg("lower_cutoff") = 0.,
        py::arg("minimum_nr_of_samples") = 0,
        py::arg("same_strand_cutoff") = 3.,
        py::arg("seed") = "",
        py::arg_v("filtering_mode",
                  EntanglementFilteringMode::OMIT_INACTIVE,
                  "EntanglementFilteringMode.OMIT_INACTIVE"),
        py::arg("n_mc_steps") = 0);

  m.def(
    "randomly_sample_entanglements",
    [](const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse,
       const pe::Box& box,
       const size_t nrOfSamples,
       const double upperCutoff,
       const double lowerCutoff,
       const size_t minimumNrOfSamples,
       const double sameStrandCutoff,
       const std::string& seed,
       const EntanglementFilteringMode filteringMode,
       const int nMcSteps) {
      return randomlyFindEntanglementsOnDemand(box,
                                               crosslinkerUniverse,
                                               nrOfSamples,
                                               upperCutoff,
                                               lowerCutoff,
                                               minimumNrOfSamples,
                                               sameStrandCutoff,
                                               seed,
                                               filteringMode,
                                               nMcSteps);
    },
    py::arg("crosslinker_universe"),
    py::arg("box"),
    py::arg("nr_of_samples"),
    py::arg("upper_cutoff"),
    py::arg("lower_cutoff") = 0.,
    py::arg("minimum_nr_of_samples") = 0,
    py::arg("same_strand_cutoff") = 3.,
    py::arg("seed") = "",
    py::arg_v("filtering_mode",
              EntanglementFilteringMode::OMIT_INACTIVE,
              "EntanglementFilteringMode.OMIT_INACTIVE"),
    py::arg("n_mc_steps") = 0);

  m.def(
    "randomly_sample_entanglements_v2",
    py::overload_cast<const pylimer_tools::entities::Universe&,
                      const size_t,
                      const double,
                      const double,
                      const size_t,
                      const double,
                      const std::string&,
                      const int,
                      const bool,
                      EntanglementFilteringMode>(&randomlyFindEntanglements),
    "Alias for randomly_sample_entanglements (v2 has been unified into "
    "standard method).",
    py::arg("universe"),
    py::arg("nr_of_samples"),
    py::arg("upper_cutoff"),
    py::arg("lower_cutoff") = 0.0,
    py::arg("minimum_nr_of_samples") = 0,
    py::arg("same_strand_cutoff") = 3.0,
    py::arg("seed") = "",
    py::arg("crosslinker_type") = 2,
    py::arg("ignore_crosslinks") = true,
    py::arg_v("filtering_mode",
              EntanglementFilteringMode::NONE,
              "EntanglementFilteringMode.NONE"));
  py::class_<pylimer_tools::topo::ResistanceDistanceAnalysis, py::smart_holder>(
    m,
    "ResistanceDistanceAnalysis",
    R"pbdoc(
      Analysis of topological resistance distances, Kirchhoff index,
      and random-walk commute times on graph Laplacians.
    )pbdoc")
    .def(py::init<const std::vector<size_t>&,
                  const std::vector<size_t>&,
                  const std::vector<double>&,
                  size_t>(),
         R"pbdoc(
           Construct from edge lists (spring_from, spring_to), optional edge weights,
           and optional number of nodes.

           :param spring_from: Source node indices of edges
           :param spring_to: Target node indices of edges
           :param weights: Optional weights (spring constants) of edges (defaults to 1.0)
           :param num_nodes: Optional total number of nodes (deduced if 0)
         )pbdoc",
         py::arg("spring_from"),
         py::arg("spring_to"),
         py::arg("weights") = std::vector<double>{},
         py::arg("num_nodes") = 0)
    .def(py::init<const Eigen::SparseMatrix<double>&>(),
         R"pbdoc(
           Construct from a pre-assembled sparse scalar graph Laplacian matrix L.

           :param L: Symmetric positive semidefinite graph Laplacian matrix
         )pbdoc",
         py::arg("L"))
    .def(py::init<const pe::Universe&>(),
         R"pbdoc(
           Construct from a Universe's atom-bond graph.

           :param universe: Universe instance
         )pbdoc",
         py::arg("universe"))
    .def_static(
      "from_matrix",
      [](const Eigen::SparseMatrix<double>& L) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::fromMatrix(L)));
      },
      R"pbdoc(
        Construct a ResistanceDistanceAnalysis instance from a sparse graph Laplacian matrix L.

        :param L: Symmetric positive semidefinite graph Laplacian matrix
        :return: ResistanceDistanceAnalysis instance
      )pbdoc",
      py::arg("L"))
    .def_static(
      "from_universe",
      [](const pe::Universe& universe) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::fromUniverse(
              universe)));
      },
      R"pbdoc(
        Construct a ResistanceDistanceAnalysis instance from a Universe's atom-bond graph.

        :param universe: Universe instance
        :return: ResistanceDistanceAnalysis instance
      )pbdoc",
      py::arg("universe"))
    .def_static(
      "create_linear_chain",
      [](size_t num_nodes) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::createLinearChain(
              num_nodes)));
      },
      R"pbdoc(
        Create a ResistanceDistanceAnalysis instance for a linear chain of num_nodes nodes.

        :param num_nodes: Number of nodes in the linear chain
        :return: ResistanceDistanceAnalysis instance
      )pbdoc",
      py::arg("num_nodes"))
    .def_static(
      "create_ring",
      [](size_t num_nodes) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::createRing(
              num_nodes)));
      },
      R"pbdoc(
        Create a ResistanceDistanceAnalysis instance for a ring graph of num_nodes nodes.

        :param num_nodes: Number of nodes in the ring (at least 3)
        :return: ResistanceDistanceAnalysis instance
      )pbdoc",
      py::arg("num_nodes"))
    .def_static(
      "create_complete",
      [](size_t num_nodes) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::createComplete(
              num_nodes)));
      },
      R"pbdoc(
        Create a ResistanceDistanceAnalysis instance for a complete graph K_N.

        :param num_nodes: Number of nodes in the complete graph (at least 2)
        :return: ResistanceDistanceAnalysis instance
      )pbdoc",
      py::arg("num_nodes"))
    .def_static(
      "create_star",
      [](size_t num_arms, size_t arm_length) {
        return std::unique_ptr<pylimer_tools::topo::ResistanceDistanceAnalysis>(
          new pylimer_tools::topo::ResistanceDistanceAnalysis(
            pylimer_tools::topo::ResistanceDistanceAnalysis::createStar(
              num_arms, arm_length)));
      },
      R"pbdoc(
        Create a ResistanceDistanceAnalysis instance for a star graph.

        :param num_arms: Number of arms (at least 2)
        :param arm_length: Number of edges per arm (defaults to 1)
        :return: ResistanceDistanceAnalysis instance
      )pbdoc",
      py::arg("num_arms"),
      py::arg("arm_length") = 1)
    .def("get_nr_of_nodes",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getNrOfNodes,
         R"pbdoc(
           Get the number of nodes in the graph.

           :return: Number of nodes
         )pbdoc")
    .def("get_nr_of_edges",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getNrOfEdges,
         R"pbdoc(
           Get the number of non-self-loop edges in the graph.

           :return: Number of edges
         )pbdoc")
    .def("get_total_edge_weight",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getTotalEdgeWeight,
         R"pbdoc(
           Get the total weight of all edges in the graph.

           :return: Total edge weight
         )pbdoc")
    .def("is_connected",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::isConnected,
         R"pbdoc(
           Check whether the graph consists of a single connected component.

           :return: True if connected, False otherwise
         )pbdoc")
    .def("get_nr_of_connected_components",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           getNrOfConnectedComponents,
         R"pbdoc(
           Get the number of connected components in the graph.

           :return: Number of connected components
         )pbdoc")
    .def(
      "get_connected_components",
      &pylimer_tools::topo::ResistanceDistanceAnalysis::getConnectedComponents,
      R"pbdoc(
        Get the list of node indices for each connected component.

        :return: List of node index lists per connected component
      )pbdoc")
    .def("get_laplacian",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getLaplacian,
         R"pbdoc(
           Get the sparse graph Laplacian matrix L.

           :return: Sparse Laplacian matrix
         )pbdoc")
    .def("get_pinned_laplacian",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::getPinnedLaplacian,
         R"pbdoc(
           Get the sparse pinned Laplacian matrix with reference nodes pinned to 1 on the diagonal.

           :return: Sparse pinned Laplacian matrix
         )pbdoc")
    .def("compute_effective_resistance",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeEffectiveResistance,
         R"pbdoc(
           Compute the effective resistance between nodes i and j.

           :param i: First node index
           :param j: Second node index
           :return: Effective resistance R_ij (or inf if in different components)
         )pbdoc",
         py::arg("i"),
         py::arg("j"))
    .def("compute_effective_resistance_batched",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeEffectiveResistanceBatched,
         R"pbdoc(
           Compute the effective resistance for a batch of node pairs using chunked multi-RHS solves.

           :param pairs: List of (i, j) node index pairs
           :param chunk_size: Number of columns per solve chunk (defaults to 128)
           :return: Vector of effective resistances for each pair
         )pbdoc",
         py::arg("pairs"),
         py::arg("chunk_size") = 128)
    .def("compute_all_pairs_effective_resistance",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeAllPairsEffectiveResistance,
         R"pbdoc(
           Compute the N x N symmetric matrix of all-pairs effective resistances.

           :param chunk_size: Number of columns per solve chunk (defaults to 128)
           :return: N x N matrix of effective resistances
         )pbdoc",
         py::arg("chunk_size") = 128)
    .def("compute_moore_penrose_pseudo_inverse",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeMoorePenrosePseudoInverse,
         R"pbdoc(
           Compute the Moore-Penrose pseudo-inverse L^+ of the graph Laplacian.

           :param chunk_size: Number of columns per solve chunk (defaults to 128)
           :return: N x N Moore-Penrose pseudo-inverse matrix
         )pbdoc",
         py::arg("chunk_size") = 128)
    .def("solve_identity_chunked",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::solveIdentityChunked,
         R"pbdoc(
           Solve L_pinned * X = I in multi-column chunks.

           :param chunk_size: Number of columns per solve chunk (defaults to 128)
           :return: N x N solution matrix X
         )pbdoc",
         py::arg("chunk_size") = 128)
    .def(
      "compute_kirchhoff_index",
      &pylimer_tools::topo::ResistanceDistanceAnalysis::computeKirchhoffIndex,
      R"pbdoc(
        Compute the Kirchhoff index K(G) = sum_{i < j} R_ij without allocating the full N x N matrix.

        :param chunk_size: Number of columns per solve chunk (defaults to 128)
        :return: Kirchhoff index (or inf if disconnected)
      )pbdoc",
      py::arg("chunk_size") = 128)
    .def("compute_kirchhoff_index_spectral",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeKirchhoffIndexSpectral,
         R"pbdoc(
           Compute the Kirchhoff index from the non-zero eigenvalues of the Laplacian.

           :return: Kirchhoff index (or inf if disconnected)
         )pbdoc")
    .def("compute_total_resistance_per_node",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeTotalResistancePerNode,
         R"pbdoc(
           Compute the sum of effective resistances from each node to all other nodes.

           :param chunk_size: Number of columns per solve chunk (defaults to 128)
           :return: Vector of total resistance per node
         )pbdoc",
         py::arg("chunk_size") = 128)
    .def("compute_resistance_centrality",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeResistanceCentrality,
         R"pbdoc(
           Compute the resistance centrality (N - 1) / sum_{j != i} R_ij for each node.

           :param chunk_size: Number of columns per solve chunk (defaults to 128)
           :return: Vector of resistance centralities
         )pbdoc",
         py::arg("chunk_size") = 128)
    .def("compute_commute_time",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::computeCommuteTime,
         R"pbdoc(
           Compute the random-walk commute time C_ij = 2 * M * R_ij between nodes i and j.

           :param i: First node index
           :param j: Second node index
           :return: Commute time C_ij
         )pbdoc",
         py::arg("i"),
         py::arg("j"))
    .def("compute_all_pairs_commute_time",
         &pylimer_tools::topo::ResistanceDistanceAnalysis::
           computeAllPairsCommuteTime,
         R"pbdoc(
           Compute the N x N matrix of all-pairs random-walk commute times.

           :param chunk_size: Number of columns per solve chunk (defaults to 128)
           :return: N x N matrix of commute times
         )pbdoc",
         py::arg("chunk_size") = 128);
};
// }

#endif
