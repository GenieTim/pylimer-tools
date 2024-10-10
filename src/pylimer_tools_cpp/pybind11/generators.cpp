#ifndef PYBIND_GENERATORS_H
#define PYBIND_GENERATORS_H

#include "../utils/MCUniverseGenerator.h"
#include "../utils/RandomWalker.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
         "Set the optimal distance between beads.",
         py::arg("distance"))
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
    .def("add_and_link_strands",
         py::overload_cast<int, std::vector<int>, double, int, double>(
           &MCUniverseGenerator::addAndLinkStrands),
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
    .def("get_universe", &MCUniverseGenerator::getUniverse, R"pbdoc(
            Fetch the current (or final) state of the universe.

            Use this method to actually retrieve the generated structure.
            )pbdoc");

  m.def("do_random_walk",
        &doRandomWalkChain,
        R"pbdoc(
            Do a random walk, return the coordinates of each point visited.
            )pbdoc",
        py::arg("chain_len"),
        py::arg("bead_distance") = 1.,
        py::arg("seed") = "");
  m.def("do_random_walk_chain_from_to",
        &doRandomWalkChainFromTo,
        R"pbdoc(
            Do a random walk from one point to another.
            )pbdoc",
        py::arg("box"),
        py::arg("from_coordinates"),
        py::arg("to_coordinates"),
        py::arg("chain_len"),
        py::arg("bead_distance") = 1.,
        py::arg("seed") = "");
  m.def("do_linear_walk_chain_from_to",
        &doLinearWalkChainFromTo,
        R"pbdoc(
            Get coordinates linearly interpolated from one point to another (both exclusive).
            )pbdoc",
        py::arg("box"),
        py::arg("from_coordinates"),
        py::arg("to_coordinates"),
        py::arg("chain_len"));
}

#endif
