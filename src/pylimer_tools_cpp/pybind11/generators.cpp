#ifndef PYBIND_GENERATORS_H
#define PYBIND_GENERATORS_H

#include "../entities/Box.h"
#include "../utils/MCUniverseGenerator.h"
#include "../utils/RandomWalker.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
using namespace pylimer_tools::utils;

void
init_pylimer_bound_generators(py::module_& m)
{

  py::class_<MCUniverseGenerator>(m, "MCUniverseGenerator", R"pbdoc(
       A :obj:`pylimer_tools_cpp.Universe` generator using a Monte-Carlo procedure.
  )pbdoc")
    .def(py::init<const double, const double, const double>(),
         py::arg("lx"),
         py::arg("ly"),
         py::arg("lz"))
    .def("set_seed",
         &MCUniverseGenerator::setSeed,
         "Set the seed for the random generator.",
         py::arg("seed"))
    .def("set_bead_distance",
         &MCUniverseGenerator::setBeadDistance,
         R"pbdoc(
         Set the mean distance between beads when doing MC stepping.
         Also used for the target cross-linker partner sampling.

         NOTE: Mainly the mean squared bead distance is effectively used in the Monte-Carlo simulation.

         :param distance: Mean distance between beads.
         :param update_mean_squared: Whether to update the mean squared distance as well, deduced from the assumed gaussian distribution in 3D (default: true).
         )pbdoc",
         py::arg("distance"),
         py::arg("update_mean_squared") = true)
    .def("get_mean_bead_distance",
         &MCUniverseGenerator::getConfiguredBeadDistance,
         "Get the currently configured mean bead distance.")
    .def("set_mean_squared_bead_distance",
         &MCUniverseGenerator::setMeanSquaredBeadDistance,
         R"pbdoc(Set the mean squared distance between beads.
         :param mean_squared_distance: Mean squared distance between beads.
         :param update_mean: Whether to update the mean bead distance as well, deduced from the assumed gaussian distribution in 3D (default: true).
         )pbdoc",
         py::arg("mean_squared_distance"),
         py::arg("update_mean") = true)
    .def("get_mean_squared_bead_distance",
         &MCUniverseGenerator::getConfiguredMeanSquaredBeadDistance,
         "Get the currently configured mean squared bead distance.")
    .def(
      "config_nr_of_mc_steps",
      &MCUniverseGenerator::configNrOfMCSteps,
      "Set the number of Monte-Carlo steps during bond length equilibration.",
      py::arg("n_steps") = 2000)
    .def("config_primary_loop_probability",
         &MCUniverseGenerator::configPrimaryLoopProbability,
         R"pbdoc(
         Configure an additional weight reduction for the primary loop formation.

         Defaults to 1., which means the general :math:`P(\vec{R})` is used without any bias.
         This results in more primary loops for shorter chains than for longer ones.

         Set to 0. to disable the formation of primary loops.
         )pbdoc",
         py::arg("probability") = 1.0)
    .def("add_crosslinkers",
         &MCUniverseGenerator::addCrosslinkers,
         R"pbdoc(
            Add the cross-linkers.

            :param nr_of_crosslinkers: Number of cross-linkers to add.
            :param crosslinker_functionality: Functionality of the cross-linkers (default: 4).
            :param crosslinker_type: Atom type of the cross-linkers (default: 2).
            :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the cross-linkers (default: true).
            )pbdoc",
         py::arg("nr_of_crosslinkers"),
         py::arg("crosslinker_functionality") = 4,
         py::arg("crosslinker_type") = 2,
         py::arg("white_noise") = true)
    .def("add_solvent_chains",
         &MCUniverseGenerator::addSolventChains,
         R"pbdoc(
            Randomly distribute additional, free chains.
            )pbdoc",
         py::arg("nr_of_solvent_chains"),
         py::arg("solvent_chain_length"),
         py::arg("solvent_atom_type") = 3,
         py::arg("white_noise") = true)
    .def("add_and_link_monofunctional_strands",
         &MCUniverseGenerator::addAndLinkMonofunctionalStrands,
         R"pbdoc(
         Add multiple monofunctional strands with specified bead types, link them to cross-links.
         )pbdoc",
         py::arg("nr_of_monofunctional_strands"),
         py::arg("monofunctional_strand_length"),
         py::arg("monofunctional_strand_atom_type") = 4)
    .def("add_and_link_strands_to_conversion",
         py::overload_cast<int, std::vector<int>, double, int, double>(
           &MCUniverseGenerator::addAndLinkStrandsToConversion),
         R"pbdoc(
            Actually add strands, link them to the previously added cross-linkers.

            :param nr_of_strands: Number of strands to add.
            :param strand_lengths: A list of integers representing the number of beads of each of the strands.
            :param crosslinker_conversion: Target conversion of cross-linkers (0: no connections to cross-links; 1: all cross-linkers fully connected).
            :param strand_atom_type: Type of atoms for the strands (default: 1).
            :param c_infinity: As needed for the end-to-end distribution, given by :math:`\langle R^2\rangle_0 = C_{\infty} N b^2`.
            )pbdoc",
         py::arg("nr_of_strands"),
         py::arg("strand_lengths"),
         py::arg("crosslinker_conversion"),
         py::arg("strand_atom_type") = 1,
         py::arg("c_infinity") = 1.)
    .def("add_and_link_strands_to_soluble_fraction",
         py::overload_cast<int, std::vector<int>, double, int, double>(
           &MCUniverseGenerator::addAndLinkStrandsToSolubleFraction),
         R"pbdoc(
            Actually add strands, link them to the previously added cross-linkers.

            :param nr_of_strands: Number of strands to add.
            :param strand_lengths: A list of integers representing the number of beads of each of the strands.
            :param soluble_fraction: Target soluble_fraction (0: no connections to cross-links; 1: all cross-linkers fully connected).
            :param strand_atom_type: Type of atoms for the strands (default: 1).
            :param c_infinity: As needed for the end-to-end distribution, given by :math:`\langle R^2\rangle_0 = C_{\infty} N b^2`.
            )pbdoc",
         py::arg("nr_of_strands"),
         py::arg("strand_lengths"),
         py::arg("soluble_fraction"),
         py::arg("strand_atom_type") = 1,
         py::arg("c_infinity") = 1.)
    .def("relax_crosslinks",
         &MCUniverseGenerator::relaxCrosslinks,
         R"pbdoc(
         Run force relaxation with the cross-linkers and their strands,
         to have the cross-links in their statistically most probable position.
         )pbdoc")
    .def("get_universe", &MCUniverseGenerator::getUniverse, R"pbdoc(
            Fetch the current (or final) state of the universe.

            Use this method to actually (MC) place beads between the cross-links and retrieve the generated structure.
            )pbdoc");

  m.def("do_random_walk",
        py::overload_cast<int, double, double, std::string>(&doRandomWalkChain),
        R"pbdoc(
            Do a random walk, return the coordinates of each point visited.
            )pbdoc",
        py::arg("chain_len"),
        py::arg("bead_distance") = 1.,
        py::arg("mean_squared_bead_distance") = 1.,
        py::arg("seed") = "");
  m.def(
    "do_random_walk_chain_from_to_mc",
    [](pe::Box& b,
       Eigen::Vector3d f,
       Eigen::Vector3d t,
       int c,
       double l,
       double l2,
       std::string s,
       int n) { return doRandomWalkChainFromToMC(b, f, t, c, l, l2, s, n); },
    R"pbdoc(
            Do a random walk from one point to another.
            Then, relax the points in between using a Metropolis-Monte Carlo simulation.
            )pbdoc",
    py::arg("box"),
    py::arg("from_coordinates"),
    py::arg("to_coordinates"),
    py::arg("chain_len"),
    py::arg("bead_distance") = 1.,
    py::arg("mean_squared_bead_distance") = 1.,
    py::arg("seed") = "",
    py::arg("n_iterations") = 1000);
  m.def(
    "do_random_walk_chain_from_to",
    [](pe::Box& b,
       Eigen::Vector3d f,
       Eigen::Vector3d t,
       int c,
       double l,
       double l2,
       std::string s) { return doRandomWalkChainFromTo(b, f, t, c, l, l2, s); },
    R"pbdoc(
            Do a random walk from one point to another.
            )pbdoc",
    py::arg("box"),
    py::arg("from_coordinates"),
    py::arg("to_coordinates"),
    py::arg("chain_len"),
    py::arg("bead_distance") = 1.,
    py::arg("mean_squared_bead_distance") = 1.,
    py::arg("seed") = "");
  m.def("do_linear_walk_chain_from_to",
        &doLinearWalkChainFromTo,
        R"pbdoc(
            Get coordinates linearly interpolated from one point to another (both exclusive).

            :param box: The box for doing PBC correction on the from/to.
            :param from_coordinates: Coordinates of the start point.
            :param to_coordinates: Coordinates of the end point.
            :param chain_len: Number of coordinates to generate between the start and end-point.
            :param include_ends: Whether to include the start and end points in the output (default: false). 
               If yes, chain_len + 2 coordinates will be returned, 
               where the first will be from_coordinates and the last will be to_coordinates.
            )pbdoc",
        py::arg("box"),
        py::arg("from_coordinates"),
        py::arg("to_coordinates"),
        py::arg("chain_len"),
        py::arg("include_ends") = false);
}

#endif
