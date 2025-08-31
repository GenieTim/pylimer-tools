#ifndef PYBIND_GENERATORS_H
#define PYBIND_GENERATORS_H

#include "../entities/Box.h"
#include "../utils/MCUniverseGenerator.h"
#include "../utils/RandomWalker.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
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
  py::class_<MaxDistanceProvider, std::shared_ptr<MaxDistanceProvider>>(
    m, "MaxDistanceProvider", R"pbdoc(
     A generic implementation of a class, that shall provide a maximum distance for the MC sampling.
     )pbdoc")
    .def("get_max_distance",
         &MaxDistanceProvider::getMaxDistance,
         R"pbdoc(
         Get the maximum distance for a given N.

         :param N: Number of segments.
         :return: Maximum distance for sampling.
         )pbdoc",
         py::arg("N"));

  py::class_<LinearMaxDistanceProvider,
             MaxDistanceProvider,
             std::shared_ptr<LinearMaxDistanceProvider>>(
    m, "LinearMaxDistanceProvider", R"pbdoc(
    For MC generation, converts the :math:`N` to a maximum distance within which to sample.
    The distance will be calculated as :math:`N \times \text{max_distance_multiplier}`.
    Useful only for performance improvements in large systems.

    :param max_distance_multiplier: Multiplier for the maximum distance.
    )pbdoc")
    .def(py::init<const double>(), py::arg("max_distance_multiplier"))
    .def("get_max_distance",
         &LinearMaxDistanceProvider::getMaxDistance,
         R"pbdoc(
         Get the maximum distance for a given N.

         :param N: Number of segments.
         )pbdoc",
         py::arg("N"));

  py::class_<ZScoreMaxDistanceProvider,
             MaxDistanceProvider,
             std::shared_ptr<ZScoreMaxDistanceProvider>>(
    m, "ZScoreMaxDistanceProvider", R"pbdoc(
     For MC generation, converts the :math:`N` to a maximum distance within which to sample.
     The distance will be calculated as :math:`\text{std_multiplier} \times \sqrt{N \times \text{in_sqrt_multiplier}}`.
     Useful only for performance improvements in large systems.

     :param std_multiplier: The Z-Score, the multiplier of the standard deviation of the end-to-end distribution.
         E.g., 3.29 for 99.9% of all conformations.
     :param in_sqrt_multiplier: The multiplier with the :math:`N` in the square root. Probably :math:`<b^2>`.
    )pbdoc")
    .def(py::init<const double, const double>(),
         py::arg("std_multiplier"),
         py::arg("in_sqrt_multiplier"))
    .def("get_max_distance",
         &ZScoreMaxDistanceProvider::getMaxDistance,
         R"pbdoc(
         Get the maximum distance for a given N.

         :param N: Number of segments.
         )pbdoc",
         py::arg("N"));

  py::class_<NoMaxDistanceProvider,
             MaxDistanceProvider,
             std::shared_ptr<NoMaxDistanceProvider>>(
    m, "NoMaxDistanceProvider", R"pbdoc(
    For MC generation, to disable the neighbour list usage.
    )pbdoc")
    .def(py::init<>())
    .def("get_max_distance",
         &NoMaxDistanceProvider::getMaxDistance,
         R"pbdoc(
         Get the maximum distance for a given N (always returns -1 to disable).

         :param N: Number of segments (ignored).
         :return: Always returns -1 to disable maximum distance checks.
         )pbdoc",
         py::arg("N"));

  py::enum_<BackTrackStatus>(m, "BackTrackStatus", R"pbdoc(
     Enum for controlling the strand linking process in linkStrandsCallback.
     )pbdoc")
    .value("STOP", BackTrackStatus::STOP, "Stop the linking process")
    .value("TRACK_FORWARD",
           BackTrackStatus::TRACK_FORWARD,
           "Continue linking forward")
    .value("TRACK_BACKWARD",
           BackTrackStatus::TRACK_BACKWARD,
           "Track backward in the linking process")
    .export_values();

  py::class_<MCUniverseGenerator>(m, "MCUniverseGenerator", R"pbdoc(
       A :obj:`pylimer_tools_cpp.Universe` generator using a Monte-Carlo procedure.

       Please cite :cite:t:`gusev_molecular_2024` and/or :cite:t:`bernhard_phantom_2025` if you use this method in your work.
  )pbdoc")
    .def(py::init<const double, const double, const double>(),
         R"pbdoc(
         Initialize a new MCUniverseGenerator with specified box dimensions.

         :param lx: Box length in x-direction (default: 10.0).
         :param ly: Box length in y-direction (default: 10.0).
         :param lz: Box length in z-direction (default: 10.0).
         )pbdoc",
         py::arg("lx") = 10.,
         py::arg("ly") = 10.,
         py::arg("lz") = 10.)
    .def("set_seed",
         &MCUniverseGenerator::setSeed,
         R"pbdoc(
         Set the seed for the random generator.

         :param seed: Random seed value for reproducible results.
         )pbdoc",
         py::arg("seed"))
    .def("set_bead_distance",
         &MCUniverseGenerator::setBeadDistance,
         R"pbdoc(
         Set the mean distance between beads when doing MC stepping.
         Also used for the target crosslinker partner sampling.

         NOTE: Mainly the mean squared bead distance is effectively used in the Monte-Carlo simulation.

         :param distance: Mean distance between beads.
         :param update_mean_squared: Whether to update the mean squared distance as well, deduced from the assumed gaussian distribution in 3D (default: true).
         )pbdoc",
         py::arg("distance"),
         py::arg("update_mean_squared") = true)
    .def("get_mean_bead_distance",
         &MCUniverseGenerator::getConfiguredBeadDistance,
         R"pbdoc(
          Get the currently configured mean bead distance.

          :return: The currently configured mean distance between beads
         )pbdoc")
    .def("set_mean_squared_bead_distance",
         &MCUniverseGenerator::setMeanSquaredBeadDistance,
         R"pbdoc(
          Set the mean squared distance between beads.

          :param mean_squared_distance: Mean squared distance between beads
          :param update_mean: Whether to update the mean bead distance as well, deduced from the assumed gaussian distribution in 3D (default: true)
         )pbdoc",
         py::arg("mean_squared_distance"),
         py::arg("update_mean") = true)
    .def("get_mean_squared_bead_distance",
         &MCUniverseGenerator::getConfiguredMeanSquaredBeadDistance,
         R"pbdoc(
          Get the currently configured mean squared bead distance.

          :return: The currently configured mean squared distance between beads
         )pbdoc")
    .def("get_nr_of_atoms",
         &MCUniverseGenerator::getCurrentNrOfAtoms,
         R"pbdoc(
          Get the current number of atoms that the universe would/will have.

          :return: Number of atoms in the generated universe
         )pbdoc")
    .def("get_nr_of_bonds",
         &MCUniverseGenerator::getCurrentNrOfBonds,
         R"pbdoc(
          Get the current number of bonds that the universe would/will have.

          :return: Number of bonds in the generated universe
         )pbdoc")
    .def("config_nr_of_mc_steps",
         &MCUniverseGenerator::configNrOfMCSteps,
         R"pbdoc(
       Set the number of Monte-Carlo steps during bond length equilibration.

       :param n_steps: Number of MC steps to perform (default: 2000)
      )pbdoc",
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
    .def("config_secondary_loop_probability",
         &MCUniverseGenerator::configSecondaryLoopProbability,
         R"pbdoc(
         Configure an additional weight reduction for the secondary loop formation.

         Defaults to 1., which means the general :math:`P(\vec{R})` is used without any bias.
         This results in more secondary loops for shorter chains than for longer ones.

         Set to 0. to disable the formation of secondary loops.
         )pbdoc",
         py::arg("probability") = 1.0)
    .def("disable_max_distance",
         &MCUniverseGenerator::disableMaxDistance,
         R"pbdoc(
         Disable the maximum distance provider to allow unlimited distance sampling.
         This may slow down performance for large systems but ensures complete sampling.
         )pbdoc")
    .def("use_linear_max_distance",
         &MCUniverseGenerator::useLinearMaxDistance,
         R"pbdoc(
          Converts the :math:`N` to a maximum distance within which to sample.
          The distance will be calculated as :math:`N \times \text{multiplier}`.
          For example, commonly, a maximum distance is :math:`N \times <b>`, 
          in which case the multiplier is :math:`<b>`.

          Useful only for performance improvements in large systems.
          )pbdoc",
         py::arg("multiplier"))
    .def("use_zscore_max_distance",
         &MCUniverseGenerator::useZScoreMaxDistance,
         R"pbdoc(
         Converts the :math:`N` to a maximum distance within which to sample.
         The distance will be calculated as :math:`\text{std_multiplier} \times \sqrt{N \times \text{in_sqrt_multiplier}}`.
         Useful only for performance improvements in large systems.

         :param std_multiplier: The Z-Score, the multiplier of the standard deviation of the end-to-end distribution.
               E.g., 3. for 99.9994% of all conformations.
         :param in_sqrt_multiplier: The multiplier with the :math:`N` in the square root. Probably :math:`<b^2>`.
         )pbdoc",
         py::arg("std_multiplier"),
         py::arg("in_sqrt_multiplier") = 1.)
    //     .def("config_max_distance_provider",
    //          &MCUniverseGenerator::configMaxDistanceProvider,
    //          R"pbdoc(
    //          For larger systems, you may want to set this value to something
    //          other than the NoMaxDistanceProvider in order to improve
    //          sampling performance.

    //          For example, you may expect a maximum end-to-end distance of a
    //          chain to be :math:`N b`, in which case you would set this value
    //          to LinearMaxDistanceProvider(:math:`b`).

    //          Note that the neighbour list is initialized when this function
    //          is called. It does respect all :math:`N` that are set at that
    //          point. Therefore, try to call this method after adding strands,
    //          before doing crosslinking. Before that, it's not relevant for
    //          performance anyway. )pbdoc", py::arg("provider"))
    .def("add_crosslinkers_at",
         &MCUniverseGenerator::addCrosslinkersAt,
         R"pbdoc(
          Add crosslinkers at specific coordinates.

          :param coordinates: Coordinates of the crosslinkers as a flat array [x1, y1, z1, x2, y2, z2, ...].
          :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
          :param crosslinker_type: Atom type of the crosslinkers (default: 2).
         )pbdoc",
         py::arg("coordinates"),
         py::arg("crosslinker_functionality") = 4,
         py::arg("crosslinker_type") = 2)
    .def("add_crosslinkers",
         &MCUniverseGenerator::addCrosslinkers,
         R"pbdoc(
            Add crosslinkers at random positions.

            :param nr_of_crosslinkers: Number of crosslinkers to add.
            :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
            :param crosslinker_type: Atom type of the crosslinkers (default: 2).
            :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).
            )pbdoc",
         py::arg("nr_of_crosslinkers"),
         py::arg("crosslinker_functionality") = 4,
         py::arg("crosslinker_type") = 2,
         py::arg("white_noise") = true)
    .def("add_randomly_functionalized_strands",
         &MCUniverseGenerator::addRandomlyFunctionalizedStrands,
         R"pbdoc(
          Add strands with randomly distributed crosslinkers in between.

          :param nr_of_strands: Number of strands to add.
          :param strand_length: Length of each strand.
          :param functionalization_probability: Proportion of beads that are made crosslink.
          :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
          :param crosslinker_type: Atom type of the crosslinkers (default: 2).
          :param strand_atom_type: Atom type of the beads that stay (default: 1).
          :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).     
     )pbdoc",
         py::arg("nr_of_strands"),
         py::arg("strand_length"),
         py::arg("functionalization_probability"),
         py::arg("crosslinker_functionality") = 4,
         py::arg("crosslinker_type") = 2,
         py::arg("strand_atom_type") = 1,
         py::arg("white_noise") = true)
    .def("add_regularly_spaced_functionalized_strands",
         &MCUniverseGenerator::addRegularlySpacedFunctionalizedStrands,
         R"pbdoc(
          Add strands with regularly spaced crosslinkers in between.

          :param nr_of_strands: Number of strands to add.
          :param strand_length: Length of each strand.
          :param spacing_between_crosslinks: Number of beads between crosslinks.
          :param offset_to_first_crosslink: Offset from start of strand to first crosslink (0-based, default: 0).
          :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
          :param crosslinker_type: Atom type of the crosslinkers (default: 2).
          :param strand_atom_type: Atom type of the beads that stay (default: 1).
          :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).     
     )pbdoc",
         py::arg("nr_of_strands"),
         py::arg("strand_length"),
         py::arg("spacing_between_crosslinks"),
         py::arg("offset_to_first_crosslink") = 0,
         py::arg("crosslinker_functionality") = 4,
         py::arg("crosslinker_type") = 2,
         py::arg("strand_atom_type") = 1,
         py::arg("white_noise") = true)
    .def("add_functionalized_strands",
         &MCUniverseGenerator::addFunctionalizedStrandsImpl,
         R"pbdoc(
          Add strands with custom crosslink placement using a selector function.

          This is a general implementation that allows for flexible crosslink placement patterns.
          The crosslink_selector function is called for each bead position and should return
          a tuple (should_convert, functionality) indicating whether that bead should be
          converted to a crosslink and what functionality it should have.

          .. tip::

             This method provides the most flexibility for creating custom crosslink patterns.
             For common patterns, consider using the more specific methods like
             :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_randomly_functionalized_strands` or
             :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_regularly_spaced_functionalized_strands`.

          Example usage::

              # Create a custom pattern with crosslinks every 5 beads starting at position 2
              def my_selector(strand_index, bead_index, total_beads):
                  if bead_index >= 2 and (bead_index - 2) % 5 == 0:
                      return (True, 4)  # Convert to crosslink with functionality 4
                  return (False, 0)     # Don't convert

              generator.add_functionalized_strands(
                  nr_of_strands=10,
                  strand_length=[20] * 10,
                  crosslink_selector=my_selector,
                  default_crosslinker_functionality=4
              )

          :param nr_of_strands: Number of strands to add.
          :param strand_length: Length of each strand (list of integers).
          :param crosslink_selector: Function that takes (strand_index, bead_index, total_beads) and returns (should_convert, functionality).
          :param default_crosslinker_functionality: Default functionality for crosslinks (used for positioning).
          :param crosslinker_type: Atom type of the crosslinkers (default: 2).
          :param strand_atom_type: Atom type of the beads that stay (default: 1).
          :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).     
     )pbdoc",
         py::arg("nr_of_strands"),
         py::arg("strand_length"),
         py::arg("crosslink_selector"),
         py::arg("default_crosslinker_functionality"),
         py::arg("crosslinker_type") = 2,
         py::arg("strand_atom_type") = 1,
         py::arg("white_noise") = true)
    .def("add_end_functionalized_strands",
         &MCUniverseGenerator::addCrosslinkStrands,
         R"pbdoc(
            Add strands with crosslinkers at the ends.

            :param nr_of_strands: Number of strands to add.
            :param strand_length: Length of each strand.
            :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
            :param crosslinker_type: Atom type of the crosslinkers (default: 2).
            :param strand_atom_type: Atom type of the beads that are not at the ends (default: 1).
            :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).
            )pbdoc",
         py::arg("nr_of_strands"),
         py::arg("strand_length"),
         py::arg("crosslinker_functionality") = 4,
         py::arg("crosslinker_type") = 2,
         py::arg("strand_atom_type") = 1,
         py::arg("white_noise") = true)
    .def("add_solvent_chains",
         &MCUniverseGenerator::addSolventChains,
         R"pbdoc(
            Randomly distribute additional, free chains.

            :param nr_of_solvent_chains: Number of solvent chains to add.
            :param solvent_chain_length: Length of each solvent chain in beads.
            :param solvent_atom_type: Atom type for solvent chain beads (default: 3).
            :param white_noise: Whether to use white noise (true) or blue noise (false) for positioning (default: true).
            )pbdoc",
         py::arg("nr_of_solvent_chains"),
         py::arg("solvent_chain_length"),
         py::arg("solvent_atom_type") = 3,
         py::arg("white_noise") = true)
    .def("add_monofunctional_strands",
         py::overload_cast<int, std::vector<int>, int>(
           &MCUniverseGenerator::addMonofunctionalStrands),
         R"pbdoc(
         Add multiple monofunctional strands with specified bead types.

         :param nr_of_monofunctional_strands: Number of monofunctional strands to add.
         :param monofunctional_strand_length: Vector specifying the length of each strand in beads.
         :param monofunctional_strand_atom_type: Atom type for the strand beads (default: 4).
         )pbdoc",
         py::arg("nr_of_monofunctional_strands"),
         py::arg("monofunctional_strand_length"),
         py::arg("monofunctional_strand_atom_type") = 4)
    .def("add_strands",
         py::overload_cast<int, std::vector<int>, int>(
           &MCUniverseGenerator::addStrands),
         R"pbdoc(
            Add strands.
            Adds them unconnected at first.
            To link them to crosslinkers, use some of the :code:`link_strand_` methods.

            :param nr_of_strands: Number of strands to add.
            :param strand_lengths: A list of integers representing the number of beads of each of the strands.
            :param strand_atom_type: Type of atoms for the strands (default: 1).
            )pbdoc",
         py::arg("nr_of_strands"),
         py::arg("strand_lengths"),
         py::arg("strand_atom_type") = 1)
    .def("add_star_crosslinkers",
         py::overload_cast<int, int, int, int, int, bool>(
           &MCUniverseGenerator::addStarCrosslinkers),
         R"pbdoc(
         Add star-like crosslinkers with pre-connected strands (useful e.g. for tetra-PEG networks).
         Each star crosslinker will have the specified functionality with strands already attached.

         .. tip::

            To have a certain polydispersity in the arms of __one__ star crosslinker,
            the stars can alternatively be created using :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_crosslinkers`,
            :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_strands` and :meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strand_to`.

         :param nr_of_stars: Number of star crosslinkers to add
         :param functionality: Functionality of each star crosslinker (number of strands)
         :param beads_per_strand: Number of beads in each strand
         :param crosslinker_atom_type: Atom type for the crosslinker
         :param strand_atom_type: Atom type for the strand beads
         :param white_noise: Whether to use white noise positioning
         )pbdoc",
         py::arg("nr_of_stars"),
         py::arg("functionality"),
         py::arg("beads_per_strand"),
         py::arg("crosslinker_atom_type") = 2,
         py::arg("strand_atom_type") = 1,
         py::arg("white_noise") = true)
    .def("link_strand",
         &MCUniverseGenerator::linkStrand,
         R"pbdoc(
          Link one particular strand to a crosslinker.
          This strand will have one bond made, to an appropriate crosslinker, 
          as chosen by the parameters associated with the strand.

          :param strand_idx: Index of the strand to be linked.
          :param c_infinity: As needed for the end-to-end distribution, given by :math:`\langle R^2\rangle_0 = C_{\infty} N b^2` (default: 1.0).
)pbdoc",
         py::arg("strand_idx"),
         py::arg("c_infinity") = 1.)
    .def("link_strand_to",
         &MCUniverseGenerator::linkStrandTo,
         R"pbdoc(
          Link a strand to a specific crosslinker.
          This assumes that you keep track of the order in which you added the crosslinkers and strands,
          as those will determine the indices.

          .. caution::

               Be aware that some few methods may change the cross-link or strand indices.

          :param strand_idx: Index of the strand to be linked.
          :param link_idx: Index of the crosslinker to be linked.
      )pbdoc",
         py::arg("strand_idx"),
         py::arg("link_idx"))
    .def("link_strands_to_conversion",
         &MCUniverseGenerator::linkStrandsToConversion,
         R"pbdoc(
            Actually link the previously added strands to the previously added crosslinkers,
            until a certain crosslink conversion is reached.

            :param crosslinker_conversion: Target conversion of crosslinkers (0: no connections to crosslinks; 1: all crosslinkers fully connected).
            :param c_infinity: As needed for the end-to-end distribution, given by :math:`\langle R^2\rangle_0 = C_{\infty} N b^2`.
            )pbdoc",
         py::arg("crosslinker_conversion"),
         py::arg("c_infinity") = 1.)
    .def("link_strands_to_soluble_fraction",
         py::overload_cast<double, double>(
           &MCUniverseGenerator::linkStrandsToSolubleFraction),
         R"pbdoc(
            Actually link the previously added strands to the previously added crosslinkers,
            until a certain soluble fraction is reached.

            :param soluble_fraction: Target soluble fraction (0: all material is soluble; 1: no material is soluble).
            :param c_infinity: As needed for the end-to-end distribution, given by :math:`\langle R^2\rangle_0 = C_{\infty} N b^2` (default: 1.0).
            )pbdoc",
         py::arg("soluble_fraction"),
         py::arg("c_infinity") = 1.)
    .def("link_strands_callback",
         &MCUniverseGenerator::linkStrandsCallback,
         R"pbdoc(
            Link strands to crosslinkers using a custom callback function to control when to stop.

            The callback function receives the current MCUniverseGenerator state and the current step number,
            and should return a BackTrackStatus value:
            - STOP: Stop the linking process
            - TRACK_FORWARD: Continue linking, make more bonds
            - TRACK_BACKWARD: Track backward in the linking process, i.e., remove bonds again

            :param linking_controller: Callback function that controls the linking process. 
                                      Function signature: (MCUniverseGenerator, int) -> BackTrackStatus
            :param c_infinity: As needed for the end-to-end distribution, given by :math:`\langle R^2\rangle_0 = C_{\infty} N b^2`.
            )pbdoc",
         py::arg("linking_controller"),
         py::arg("c_infinity") = 1.)
    .def(
      "link_strand_to_strand",
      py::overload_cast<const double>(&MCUniverseGenerator::linkStrandToStrand),
      R"pbdoc(
         Link two free strand ends together directly.
         This method finds free strand ends and combines them into a single strand based on
         end-to-end distance probability distributions.

         .. caution::

               - The strand indices will change after this method is called.
               - The probability of linking to the end of a free strand is the same as 
                 linking to an end with a distance 0.

         :param c_infinity: Statistical parameter for end-to-end distance calculation
         :return: True if a successful link was made, False if no suitable pair was found
         )pbdoc",
      py::arg("c_infinity") = 1.0)
    .def("link_strands_to_strands_to_conversion",
         &MCUniverseGenerator::linkStrandsToStrandsToConversion,
         R"pbdoc(
         Link free strand ends to each other until target conversion is reached.
         This method repeatedly calls :meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strand_to_strand` 
         in a slightly more efficient manner than you could from Python,
         until the specified conversion of free strand ends is achieved.

         .. warning::

            The strands are merged during this process, i.e., from two strands results one strand.
            Consequently, calling :meth:`~pylimer_tools_cpp.MCUniverseGenerator.get_current_strands_conversion` after this method will not return
            the requested `target_strand_conversion`, since that one is taken relative to the number of strands when calling this method.

         .. caution::

            Beware of the notes on :meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strand_to_strand`.

         :param target_strand_conversion: Target conversion of free strand ends (0.0 to 1.0)
         :param c_infinity: Statistical parameter for end-to-end distance calculation
         )pbdoc",
         py::arg("target_strand_conversion"),
         py::arg("c_infinity") = 1.0)
    .def("remove_soluble_fraction",
         &MCUniverseGenerator::removeSolubleFraction,
         R"pbdoc(
            Remove soluble fraction (as determined by phantom force relaxation) of the strands and crosslinks.

            :param rescale_box: Whether to rescale the box dimensions to maintain constant density (default: true).
            )pbdoc",
         py::arg("rescale_box") = true)
    .def("relax_crosslinks",
         &MCUniverseGenerator::relaxCrosslinks,
         R"pbdoc(
         Run force relaxation with the crosslinkers and their strands,
         to have the crosslinks in their statistically most probable position.
         )pbdoc")
    .def("get_current_nr_of_atoms",
         &MCUniverseGenerator::getCurrentNrOfAtoms,
         R"pbdoc(
         Get the current number of atoms that the universe would/will have.

         :return: Number of atoms in the generated universe.
         )pbdoc")
    .def("get_current_nr_of_bonds",
         &MCUniverseGenerator::getCurrentNrOfBonds,
         R"pbdoc(
         Get the current number of bonds that the universe would/will have.

         :return: Number of bonds in the generated universe.
         )pbdoc")
    .def("get_current_crosslinker_conversion",
         &MCUniverseGenerator::getCurrentCrosslinkerConversion,
         R"pbdoc(
         Get the current conversion of crosslinkers, i.e., the fraction of crosslinkers
         that have been linked to strands.

         :return: Crosslinker conversion as a fraction between 0 and 1.
    )pbdoc")
    .def("get_current_strand_conversion",
         &MCUniverseGenerator::getCurrentStrandsConversion,
         R"pbdoc(
         Get the current conversion of strands, i.e., the fraction of strand ends that have been consumed.

         :return: Strand conversion as a fraction between 0 and 1.
         )pbdoc")
    .def("validate",
         &MCUniverseGenerator::validateInternalState,
         R"pbdoc(
         Check whether the internal state of the generator is valid.
         Throws errors if not. 
         Should in principle always be valid when called from Python; if not, there is a bug in the code.
         )pbdoc")
    .def("get_universe",
         &MCUniverseGenerator::getUniverse,
         R"pbdoc(
            Fetch the current (or final) state of the universe.

            Use this method to actually (MC) place beads between the crosslinks and retrieve the generated structure.

            :return: The generated Universe object containing all atoms, bonds, and their coordinates.
            )pbdoc")
    .def("get_force_relaxation",
         &MCUniverseGenerator::getForceRelaxation,
         R"pbdoc(
         Get an instance of the force relaxation procedure.
         This is a useful shorthand e.g. to skip the sampling of beads within in the strands.

         :return: Configured MEHPForceRelaxation instance.
         )pbdoc")
    .def("get_force_balance",
         &MCUniverseGenerator::getForceBalance,
         R"pbdoc(
         Get an instance of the force balance procedure.
         This is a useful shorthand e.g. to skip the sampling of beads within in the strands.

         :return: Configured MEHPForceBalance instance.
         )pbdoc")
    .def("get_force_balance2",
         &MCUniverseGenerator::getForceBalance2,
         R"pbdoc(
              Get an instance of the force balance procedure.
              This is a useful shorthand e.g. to skip the sampling of beads within in the strands.

              :return: Configured MEHPForceBalance2 instance.
              )pbdoc")
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const MCUniverseGenerator& gen) {
        return py::make_tuple(pylimer_tools::utils::serializeToString(gen));
      },
      [](py::tuple t) {
        std::string in = t[0].cast<std::string>();
        MCUniverseGenerator gen;
        pylimer_tools::utils::deserializeFromString(gen, in);
        return gen;
      }))
#endif
    .def("__copy__", [](const MCUniverseGenerator& gen) {
      return MCUniverseGenerator(gen);
    });

  m.def("do_random_walk",
        py::overload_cast<int, double, double, std::string>(&doRandomWalkChain),
        R"pbdoc(
            Do a random walk, return the coordinates of each point visited.

            :param chain_len: Length of the chain to generate.
            :param bead_distance: Mean distance between consecutive beads (default: 1.0).
            :param mean_squared_bead_distance: Mean squared distance between consecutive beads (default: 1.0).
            :param seed: Random seed for reproducibility (default: empty string for random seed).
            :return: Coordinates of each point as a flat array (i.e., [x1, y1, z1, x2, y2, z2, ...]).
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
       const int c,
       const double l,
       const double l2,
       std::string s,
       const int n) {
      return doRandomWalkChainFromToMC(b, f, t, c, l, l2, s, n);
    },
    R"pbdoc(
            Do a random walk from one point to another.
            Then, relax the points in between using a Metropolis-Monte Carlo simulation.

            :param box: Simulation box for periodic boundary conditions.
            :param from_coordinates: Starting coordinates as 3D vector.
            :param to_coordinates: Target coordinates as 3D vector.
            :param chain_len: Number of beads to place between start and end points.
            :param bead_distance: Mean distance between consecutive beads (default: 1.0).
            :param mean_squared_bead_distance: Mean squared distance between consecutive beads (default: 1.0).
            :param seed: Random seed for reproducibility (default: empty string for no seed).
            :param n_iterations: Number of Monte Carlo iterations for relaxation (default: 10000).
            :return: Coordinates of the relaxed chain as a flat array (i.e., [x1, y1, z1, x2, y2, z2, ...]).
            )pbdoc",
    py::arg("box"),
    py::arg("from_coordinates"),
    py::arg("to_coordinates"),
    py::arg("chain_len"),
    py::arg("bead_distance") = 1.,
    py::arg("mean_squared_bead_distance") = 1.,
    py::arg("seed") = "",
    py::arg("n_iterations") = 10000);
  m.def(
    "do_random_walk_chain_from_to",
    [](pe::Box& b,
       Eigen::Vector3d f,
       Eigen::Vector3d t,
       const int c,
       const double l,
       const double l2,
       std::string s) { return doRandomWalkChainFromTo(b, f, t, c, l, l2, s); },
    R"pbdoc(
            Do a random walk from one point to another.

            :param box: Simulation box for periodic boundary conditions.
            :param from_coordinates: Starting coordinates as 3D vector.
            :param to_coordinates: Target coordinates as 3D vector.
            :param chain_len: Number of beads to place between start and end points.
            :param bead_distance: Mean distance between consecutive beads (default: 1.0).
            :param mean_squared_bead_distance: Mean squared distance between consecutive beads (default: 1.0).
            :param seed: Random seed for reproducibility (default: empty string for no seed).
            :return: Coordinates of the chain as a flat array (i.e., [x1, y1, z1, x2, y2, z2, ...]).
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
