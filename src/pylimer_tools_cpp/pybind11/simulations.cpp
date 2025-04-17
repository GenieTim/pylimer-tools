#ifndef PYBIND_SIM_H
#define PYBIND_SIM_H

#include "../entities/Universe.h"
#include "../sim/DPDSimulator.h"
#include "../sim/MEHPForceBalance.h"
#include "../sim/MEHPForceBalance2.h"
#include "../sim/MEHPForceEvaluator.h"
#include "../sim/MEHPForceRelaxation.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;

using namespace pylimer_tools::sim;

namespace pylimer_tools::sim::mehp {
class PyMEHPForceEvaluator : public MEHPForceEvaluator
{
public:
  using MEHPForceEvaluator::getNetwork;

  /* Trampoline */
  virtual double evaluateStressContribution(double springDistances[3],
                                            size_t i,
                                            size_t j,
                                            size_t springIndex) const override
  {
    PYBIND11_OVERRIDE_PURE(double,
                           /* Return type */
                           MEHPForceEvaluator,
                           /* Parent class */
                           evaluateStressContribution,
                           /* Name of function in C++ */
                           springDistances,
                           i,
                           j,
                           springIndex /* Arguments */
    );
  }

  typedef std::pair<double, std::vector<double>> returntype;
  /* Trampoline */
  virtual returntype evaluateForceAndGradient(
    const size_t n,
    const Eigen::VectorXd& springDistances,
    bool requiresGradient) const
  {
    PYBIND11_OVERRIDE_PURE(
      returntype,
      /* Return type */
      MEHPForceEvaluator,
      /* Parent class */
      evaluateForceSetGradient,
      /* Name of function in C++
                                                       (must
         match Python name) */
      n,
      springDistances,
      requiresGradient /* Argument(s) */
    );
  }

  // actually overriding function, but simplifying for python possibilities
  double evaluateForceSetGradient(const size_t n,
                                  const Eigen::VectorXd& springDistances,
                                  double* grad) const override
  {
    std::pair<double, std::vector<double>> trampolineResult =
      this->evaluateForceAndGradient(n, springDistances, grad != nullptr);
    if (grad != nullptr) {
      assert(trampolineResult.second.size() == n);
      for (size_t i = 0; i < n; ++i) {
        grad[i] = trampolineResult.second[i];
      }
    }
    return trampolineResult.first;
  }

  void prepareForEvaluations() override {};
};
}

void
init_pylimer_bound_sim(py::module_& m)
{
  ////////////////////////////////////////////////////////////////
  // MARK: Output Quantities

  py::enum_<ComputedDoubleValues>(m, "ComputedDoubleValues")
    .value("TIMESTEP", ComputedDoubleValues::TIMESTEP)
    .value("TIME", ComputedDoubleValues::TIME)
    .value("VOLUME", ComputedDoubleValues::VOLUME)
    .value("PRESSURE", ComputedDoubleValues::PRESSURE)
    .value("TEMPERATURE", ComputedDoubleValues::TEMPERATURE)
    .value("STRESS_XX", ComputedDoubleValues::STRESS_XX)
    .value("STRESS_YY", ComputedDoubleValues::STRESS_YY)
    .value("STRESS_ZZ", ComputedDoubleValues::STRESS_ZZ)
    .value("STRESS_XY", ComputedDoubleValues::STRESS_XY)
    .value("STRESS_XZ", ComputedDoubleValues::STRESS_XZ)
    .value("STRESS_YZ", ComputedDoubleValues::STRESS_YZ)
    .value("STRESS_NXY", ComputedDoubleValues::STRESS_NXY)
    .value("STRESS_NXZ", ComputedDoubleValues::STRESS_NXZ)
    .value("STRESS_NYZ", ComputedDoubleValues::STRESS_NYZ)
    .value("GAMMA", ComputedDoubleValues::GAMMA)
    .value("RESIDUAL", ComputedDoubleValues::RESIDUAL)
    .value("MEAN_B", ComputedDoubleValues::MEAN_B)
    .value("MAX_B", ComputedDoubleValues::MAX_B)
    .value("MSD", ComputedDoubleValues::MSD);

  py::enum_<ComputedIntValues>(m, "ComputedIntValues")
    .value("STEP", ComputedIntValues::STEP)
    .value("NUM_SHIFT", ComputedIntValues::NUM_SHIFT)
    .value("NUM_RELOC", ComputedIntValues::NUM_RELOC)
    .value("NUM_ATOMS", ComputedIntValues::NUM_ATOMS)
    .value("NUM_EXTRA_ATOMS", ComputedIntValues::NUM_EXTRA_ATOMS)
    .value("NUM_BONDS", ComputedIntValues::NUM_BONDS)
    .value("NUM_EXTRA_BONDS", ComputedIntValues::NUM_EXTRA_BONDS)
    .value("NUM_BONDS_TO_FORM", ComputedIntValues::NUM_BONDS_TO_FORM);

  py::class_<OutputConfiguration>(m, "OutputConfiguration", py::module_local())
    .def(py::init<>(), "Get an instance of this struct")
    .def_readwrite("int_values", &OutputConfiguration::intValues)
    .def_readwrite("double_values", &OutputConfiguration::doubleValues)
    .def_readwrite("use_every",
                   &OutputConfiguration::useEvery,
                   R"pbdoc(
     For autocorrelation/averaging, how often to include values
     )pbdoc")
    .def_readwrite("append",
                   &OutputConfiguration::append,
                   R"pbdoc(
     Whether to append to the file or truncate it
     )pbdoc")
    .def_readwrite("filename",
                   &OutputConfiguration::filename,
                   R"pbdoc(
     The file to write to. Empty means standard output (console).
     )pbdoc")
    .def_readwrite("output_every",
                   &OutputConfiguration::outputEvery,
                   R"pbdoc(
     How often to write the values to the output.
     For averages, this value also says how many values will be averaged.
     )pbdoc");

  /**
   * ////////////////////////////////////////////////////////////////
   * MEHP
   * ////////////////////////////////////////////////////////////////
   */

  py::enum_<mehp::ExitReason>(m, "ExitReason")
    .value("UNSET", mehp::ExitReason::UNSET)
    .value("MAX_STEPS", mehp::ExitReason::MAX_STEPS)
    .value("F_TOLERANCE", mehp::ExitReason::F_TOLERANCE)
    .value("X_TOLERANCE", mehp::ExitReason::X_TOLERANCE)
    .value("FAILURE", mehp::ExitReason::FAILURE)
    .value("OTHER", mehp::ExitReason::OTHER);

  m.def("inverse_langevin",
        &mehp::langevin_inv,
        R"pbdoc(
     A somewhat accurate (for :math:`x \in (-1, 1)`) implementation of the inverse Langevin.

     Source: https://scicomp.stackexchange.com/a/30251
  )pbdoc",
        py::arg("x"));

  ////////////////////////////////////////////////////////////////
  // MARK: Network structures
  py::class_<mehp::Network>(m,
                            "SimplifiedNetwork",
                            R"pbdoc(
     A more efficient structure of the network for use in MEHP.
     Consists usually only of the crosslinkers.
 )pbdoc")
    .def_readonly("box_lengths", &mehp::Network::L)
    .def_readonly("volume", &mehp::Network::vol)
    .def_readonly("nr_of_nodes", &mehp::Network::nrOfNodes)
    .def_readonly("nr_of_crosslinks", &mehp::Network::nrOfNodes)
    .def_readonly("nr_of_springs", &mehp::Network::nrOfSprings)
    // .def_readonly("nrOfLoops", &mehp::Network::nrOfLoops)
    .def_readonly("coordinates", &mehp::Network::coordinates)
    .def_readonly("old_atom_ids", &mehp::Network::oldAtomIds)
    .def_readonly("spring_coordinate_index_a",
                  &mehp::Network::springCoordinateIndexA)
    .def_readonly("spring_coordinate_index_b",
                  &mehp::Network::springCoordinateIndexB)
    .def_readonly("spring_index_a", &mehp::Network::springIndexA)
    .def_readonly("spring_index_b", &mehp::Network::springIndexB)
    // .def_readonly("springIsActive", &mehp::Network::springIsActive)
    ;

  py::class_<mehp::ForceBalanceNetwork>(m,
                                        "SimplifiedBalanceNetwork",
                                        R"pbdoc(
     A more efficient structure of the network for use in MEHP force balance.
     Consists usually only of the cross- and slip-links.

     Assumed terminology: a spring is approximately a strand/chain,
     whereas a partial spring is the spring between a cross-link and an entanglement-link (slip-link).
 )pbdoc")
    .def_readonly("box_lengths", &mehp::ForceBalanceNetwork::L)
    .def_readonly("volume", &mehp::ForceBalanceNetwork::vol)
    .def_readonly("nr_of_crosslinks", &mehp::ForceBalanceNetwork::nrOfNodes)
    .def_readonly("nr_of_links", &mehp::ForceBalanceNetwork::nrOfLinks)
    .def_readonly("nr_of_springs", &mehp::ForceBalanceNetwork::nrOfSprings)
    .def_readonly("nr_of_partial_springs",
                  &mehp::ForceBalanceNetwork::nrOfPartialSprings)
    // .def_readonly("nrOfLoops", &mehp::Network::nrOfLoops)
    .def_readonly("coordinates", &mehp::ForceBalanceNetwork::coordinates)
    .def_readonly("old_atom_ids", &mehp::ForceBalanceNetwork::oldAtomIds)
    .def_readonly("spring_coordinate_index_a",
                  &mehp::ForceBalanceNetwork::springCoordinateIndexA)
    .def_readonly("spring_coordinate_index_b",
                  &mehp::ForceBalanceNetwork::springCoordinateIndexB)
    .def_readonly("spring_part_coordinate_index_a",
                  &mehp::ForceBalanceNetwork::springPartCoordinateIndexA)
    .def_readonly("spring_part_coordinate_index_b",
                  &mehp::ForceBalanceNetwork::springPartCoordinateIndexB)
    .def_readonly("spring_index_a", &mehp::ForceBalanceNetwork::springIndexA)
    .def_readonly("spring_index_b", &mehp::ForceBalanceNetwork::springIndexB)
    .def_readonly("spring_part_index_a",
                  &mehp::ForceBalanceNetwork::springPartIndexA)
    .def_readonly("spring_part_index_b",
                  &mehp::ForceBalanceNetwork::springPartIndexB)
    .def_readonly("link_is_sliplink",
                  &mehp::ForceBalanceNetwork::linkIsSliplink)
    .def_readonly("local_to_global_spring_index",
                  &mehp::ForceBalanceNetwork::localToGlobalSpringIndex)
    .def_readonly("spring_indices_of_links",
                  &mehp::ForceBalanceNetwork::springIndicesOfLinks)
    .def_readonly("link_indices_of_springs",
                  &mehp::ForceBalanceNetwork::linkIndicesOfSprings)
    .def_readonly("nr_of_crosslink_swaps_endured",
                  &mehp::ForceBalanceNetwork::nrOfCrosslinkSwapsEndured)
    .def_readonly("spring_contour_length",
                  &mehp::ForceBalanceNetwork::springsContourLength)
    .def_readonly("partial_to_full_spring_index",
                  &mehp::ForceBalanceNetwork::partialToFullSpringIndex)
    .def_readonly("spring_part_box_offset",
                  &mehp::ForceBalance2Network::springPartBoxOffset)
    // .def_readonly("springIsActive", &mehp::Network::springIsActive)
    ;

  py::class_<mehp::ForceBalance2Network>(m,
                                         "SimplifiedBalance2Network",
                                         R"pbdoc(
A more efficient structure of the network for use in MEHP force balance 2.
Consists usually only of the cross- and entanglement-links.

The terminology is a bit more consistent here:
the strands are the chains between two junctions ("nodes"), whereas the 
links will be both these nodes as well as entanglement-links, 
and finally springs will be any number of connected bonds between links.
)pbdoc")
    .def_readonly("box_lengths", &mehp::ForceBalance2Network::L)
    .def_readonly("nr_of_crosslinks", &mehp::ForceBalance2Network::nrOfNodes)
    .def_readonly("nr_of_links", &mehp::ForceBalance2Network::nrOfLinks)
    .def_readonly("nr_of_strands", &mehp::ForceBalance2Network::nrOfStrands)
    .def_readonly("nr_of_springs", &mehp::ForceBalance2Network::nrOfSprings)
    // .def_readonly("nrOfLoops", &mehp::Network::nrOfLoops)
    .def_readonly("coordinates", &mehp::ForceBalance2Network::coordinates)
    .def_readonly("old_atom_ids", &mehp::ForceBalance2Network::oldAtomIds)
    .def_readonly("old_atom_types", &mehp::ForceBalance2Network::oldAtomTypes)
    .def_readonly("spring_coordinate_index_a",
                  &mehp::ForceBalance2Network::springCoordinateIndexA)
    .def_readonly("spring_coordinate_index_b",
                  &mehp::ForceBalance2Network::springCoordinateIndexB)
    .def_readonly("spring_index_a", &mehp::ForceBalance2Network::springIndexA)
    .def_readonly("spring_index_b", &mehp::ForceBalance2Network::springIndexB)
    .def_readonly("link_is_entanglement",
                  &mehp::ForceBalance2Network::linkIsEntanglement)
    .def_readonly("spring_contour_length",
                  &mehp::ForceBalance2Network::springsContourLength)
    .def_readonly("spring_indices_of_strand",
                  &mehp::ForceBalance2Network::springIndicesOfStrand)
    .def_readonly("strand_indices_of_link",
                  &mehp::ForceBalance2Network::strandIndicesOfLink)
    .def_readonly("link_indices_of_strand",
                  &mehp::ForceBalance2Network::linkIndicesOfStrand)
    .def_readonly("nr_of_crosslink_swaps_endured",
                  &mehp::ForceBalance2Network::nrOfCrosslinkSwapsEndured)
    .def_readonly("strand_index_of_spring",
                  &mehp::ForceBalance2Network::strandIndexOfSpring)
    .def_readonly("spring_box_offset",
                  &mehp::ForceBalance2Network::springBoxOffset);

  ////////////////////////////////////////////////////////////////
  // MARK: Force evaluators
  py::class_<mehp::MEHPForceEvaluator, mehp::PyMEHPForceEvaluator>(
    m,
    "MEHPForceEvaluator",
    R"pbdoc(
     The base interface to change the way the force is evaluated during a MEHP run.
    )pbdoc")
    .def(py::init<>())
    .def_property_readonly("network", &mehp::MEHPForceEvaluator::getNetwork)
    .def_property("is_2d",
                  &mehp::MEHPForceEvaluator::getIs2D,
                  &mehp::MEHPForceEvaluator::setIs2D)
    //     .def("evaluateForceSetGradient",
    //          py::overload_cast<const size_t,
    //                            const Eigen::VectorXd&,
    //                            const Eigen::VectorXd&,
    //                            bool>(
    //            &mehp::MEHPForceEvaluator::evaluateForceSetGradient))
    .def("evaluate_stress_contribution",
         &mehp::MEHPForceEvaluator::evaluateStressContribution,
         R"pbdoc(
          An evaluation of the stress-contribution.

          :param springDistances: the three coordinate differences for one spring.
          :param i: the row index of the stress tensor
          :param j: the column index of the stress tensor
    )pbdoc",
         py::arg("spring_distances"),
         py::arg("i"),
         py::arg("j"),
         py::arg("spring_index"));

  //   py::class_<mehp::PyMEHPForceEvaluator, mehp::MEHPForceEvaluator>(
  //     m, "CustomMEHPForceEvaluator", R"pbdoc(
  //      The Python access to implement a custom force to be evaluated during a
  //      MEHP run.
  //     )pbdoc")
  //     .def(py::init<>())
  //     .def("evaluateForceAndGradient",
  //          &mehp::PyMEHPForceEvaluator::evaluateForceAndGradient,
  //          R"pbdoc(
  //      One of the two functions to override, the other being
  //      :func:`~pylimer_tools_cpp.MEHPForceEvaluator.evaluateStressContribution`.

  //      :param n: the dimensionality of the problem (the nr. of spring
  //      coordinates) :param springDistances: the sequential (x, y, z) spring
  //      distances :param displacements: the displacements from the original
  //      coordinates
  //           (accessible by
  //           :func:`~pylimer_tools_cpp.CustomMEHPForceEvaluator.getNetwork().coordinates`)
  //      :param gradientNeeded: whether the gradient should be computed and
  //      returned

  //      Returns:
  //           - force: the result of the force computation.
  //           - gradient: the result of the gradient computation.
  //                Only needed if the parameter `gradientNeeded` is true,
  //                otherwise an empty list is sufficient.
  //     )pbdoc",
  //          py::arg("n"),
  //          py::arg("springDistances"),
  //          py::arg("displacements"),
  //          py::arg("gradientNeeded"));

  py::class_<mehp::SimpleSpringMEHPForceEvaluator, mehp::MEHPForceEvaluator>(
    m,
    "SimpleSpringMEHPForceEvaluator",
    R"pbdoc(
     This is equal to a spring evaluator for Gaussian chains.

     The force for a certain spring is given by:
     :math:`f = 0.5 \cdot \kappa r`,
     where :math:`r` is the spring [between cross-linkers] length.

     Recommended optimization algorithm: "LD_LBFGS"

     :param kappa: the spring constant :math:`\kappa`
    )pbdoc")
    .def(py::init<double>(), py::arg("kappa") = 1.0);

  py::class_<mehp::NonGaussianSpringForceEvaluator, mehp::MEHPForceEvaluator>(
    m,
    "NonGaussianSpringForceEvaluator",
    R"pbdoc(
     This is equal to a spring evaluator for Langevin chains.

     The force for a certain spring is given by:
     :math:`f = 0.5 \cdot \frac{1}{l} \scriptL^{-1}(\frac{r}{N\cdot l})`,
     where :math:`r` is the spring [between cross-linkers] length
     and :math:`\scriptL^{-1}` the inverse langevin function.

     Please note that the inverse langevin is only approximated.

     Recommended optimization algorithm: "LD_MMA"

     :param kappa: the spring constant :math:`\kappa`
     :param N: The number of links in a spring
     :param l: The  the length of a spring in the chain
    )pbdoc")
    .def(py::init<double, double, double>(),
         "Initialize this ForceEvaluator",
         py::arg("kappa") = 1.0,
         py::arg("N") = 1.0,
         py::arg("l") = 1.0);

  ////////////////////////////////////////////////////////////////
  // MARK: Force Relaxation
  py::class_<mehp::MEHPForceRelaxation>(m,
                                        "MEHPForceRelaxation",
                                        R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.

    This is the first of three force relaxation methods available in this library.
    The relevant feature of this implementation is the configurable spring potential.
    Consequently, it offers a variety of configurable non-linear solvers using NLoptLib.
     )pbdoc")
    .def(py::init<pe::Universe,
                  int,
                  bool,
                  mehp::MEHPForceEvaluator*,
                  double,
                  bool,
                  bool>(),
         R"pbdoc(
          Instantiate the simulator for a certain universe.

          :param universe: the universe to simulate with
          :param crosslinker_type: The atom type of the cross-linkers. Needed to reduce the network.
          :param is2d: Whether to ignore the z direction.
          :param force_evaluator: The force evaluator to use
          :param kappa: The spring constant
          :param remove_2functional_crosslinkers: Whether to replace two-functional crosslinkers with a "normal" chain bead
          :param remove_dangling_chains: Whether to remove dangling chains before running the simulation.
               **Caution*: Removing the dangling chains will result in incorrect results fo the computation of
               :func:`~pylimer_tools_cpp.MEHPForceRelaxation.getSolubleWeightFraction()` and
               :func:`~pylimer_tools_cpp.MEHPForceRelaxation.getDanglingWeightFraction()`
          )pbdoc",
         py::arg("universe"),
         py::arg("crosslinker_type") = 2,
         py::arg("is_2d") = false,
         py::arg("force_evaluator") = nullptr,
         py::arg("kappa") = 1.0,
         py::arg("remove_2functional_crosslinkers") = false,
         py::arg("remove_dangling_chains") = false)
    .def("run_force_relaxation",
         &mehp::MEHPForceRelaxation::runForceRelaxation,
         R"pbdoc(
          Run the simulation.
          Note that the final state of the minimization is persisted and reused if you use this method again.
          This is useful if you want to run a global optimization first and add a local one afterwards.
          As a consequence though, you cannot simply benchmark only this method; you must include the setup.

          :param algorithm: The algorithm to use for the force relaxation. Choices: see `NLopt Algorithms <https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/>`_
          :param maxNrOfSteps: The maximum number of steps to do during the simulation.
          :param xTolerance: The tolerance of the displacements as an exit condition.
          :param fTolerance: The tolerance of the force as an exit condition.
          :param is2d: Specify true if you want to evaluate the force relation only in x and y direction.
          )pbdoc",
         py::arg("algorithm") = "LD_MMA",
         py::arg("max_nr_of_steps") = 250000,
         py::arg("x_tolerance") = 1e-12,
         py::arg("f_tolerance") = 1e-9)
    // .def("getForceEvaluator", &mehp::MEHPForceRelaxation::getForceEvaluator,
    // R"pbdoc(
    //      Query the currently used force evaluator.
    // )pbdoc")
    .def("set_force_evaluator",
         &mehp::MEHPForceRelaxation::setForceEvaluator,
         R"pbdoc(
          Reset the currently used force evaluator.
     )pbdoc",
         py::arg("force_evaluator"))
    .def("config_rerun_epsilon",
         &mehp::MEHPForceRelaxation::configRerunEps,
         R"pbdoc(
          Configure the offset from the lower and upper bounds for the simulation to suggest another run (
               See: :func:`~pylimer_tools_cpp.MEHPForceRelaxation.requiresAnotherRun()`
          ).
         )pbdoc",
         py::arg("epsilon") = 1e-3)
    .def("config_step_output",
         &mehp::MEHPForceRelaxation::configStepOutput,
         R"pbdoc(
          Set which values to log.

          Arguments:
               - values: a list of OutputConfiguration structs
     )pbdoc",
         py::arg("output_configuration"))
    .def("assume_box_large_enough",
         &mehp::MEHPForceRelaxation::configAssumeBoxLargeEnough,
         R"pbdoc(
          Configure whether to run PBC on the bonds or not.

          If your bonds could get larger than half the box length, this must be kept false (default).
          Otherwise, you can set it to true and therewith get some securities.
         )pbdoc",
         py::arg("box_large_enough") = false)
    .def("get_force",
         &mehp::MEHPForceRelaxation::getForce,
         R"pbdoc(
          Returns the force at the current state of the simulation.
     )pbdoc")
    .def("get_residuals",
         &mehp::MEHPForceRelaxation::getResiduals,
         R"pbdoc(
          Returns the residuals at the current state of the simulation.
     )pbdoc")
    .def("get_residual_norm",
         &mehp::MEHPForceRelaxation::getResidualNorm,
         R"pbdoc(
          Returns the residual norm at the current state of the simulation.
     )pbdoc")
    .def("get_pressure",
         &mehp::MEHPForceRelaxation::getPressure,
         R"pbdoc(
          Returns the pressure at the current state of the simulation.
     )pbdoc")
    .def("get_stress_tensor",
         &mehp::MEHPForceRelaxation::getStressTensor,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc")
    .def("get_gamma_factors",
         &mehp::MEHPForceRelaxation::getGammaFactors,
         R"pbdoc(
          Computes the gamma factor for each spring as part of the ANT/MEHP formulism.

          :math:`\gamma_{\eta} = \frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`, with (here)
          :math:`R_{0,\eta}^2 = N_\eta \cdot ` the parameter `b0_squared`.
          You can obtain this parameter e.g. by doing melt simulations at different lengths,
          it's the slope you obtain.

          :param b0_squared: Part of the denominator in the equation of :math:`\Gamma`.
               If :math:`-1.0` (default), the network is used for determination (which is not accurate), the system is assumed to be phantom.
               For real systems, the value could be determined by :func:`~pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance()`
               on the melt system, with subsequent division by the nr of bonds in the chain.

          See also :func:`~pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor` for the mean of these.
         )pbdoc",
         py::arg("b0_squared") = -1.0)
    .def("get_gamma_factor",
         &mehp::MEHPForceRelaxation::getGammaFactor,
         R"pbdoc(
          Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

          :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
          which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
          where :math:`\eta` is the index of a particular strand,
          :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`= N_{\eta} b^2$`,
          :math:`N_{\eta}` is the number of atoms in this strand :math:`\eta`,
          :math:`b` its mean square bond length,
          :math:`\nu` the volume fraction,
          :math:`T` the temperature and
          :math:`k_B` Boltzmann's constant.

          :param b0_squared: Part of the denominator in the equation of :math:`\Gamma`.
               If :math:`-1.0` (default), the network is used for determination (which is not accurate), the system is assumed to be phantom.
               For real systems, the value could be determined by :func:`~pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance()`
               on the melt system, with subsequent division by the nr of bonds in the chain.
          :param nr_of_chains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains.
     )pbdoc",
         py::arg("b0_squared") = -1.0,
         py::arg("nr_of_chains") = -1)
    .def("get_nr_of_nodes",
         &mehp::MEHPForceRelaxation::getNrOfNodes,
         R"pbdoc(
           Get the number of nodes considered in this simulation.
     )pbdoc")
    .def("get_nr_of_springs",
         &mehp::MEHPForceRelaxation::getNrOfSprings,
         R"pbdoc(
          Get the number of springs considered in this simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc")
    .def("get_ids_of_active_nodes",
         &mehp::MEHPForceRelaxation::getIdsOfActiveNodes,
         R"pbdoc(
          Get the atom ids of the nodes that are considered active.

          Arguments:
           - :param tolerance: springs under this length are considered inactive. A node is active if it has > 2 active springs.
     )pbdoc",
         py::arg("tolerance") = 1e-3,
         py::arg("minimum_nr_of_active_connections") = 2,
         py::arg("maximum_nr_of_active_connections") = -1)
    .def("get_nr_of_active_nodes",
         &mehp::MEHPForceRelaxation::getNrOfActiveNodes,
         R"pbdoc(
           Get the number of active nodes remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive.
          :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
          :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
               Use a value < 0 to indicate that there is no maximum number of active connections.
          :param usePartial: Whether to use the partial spring distances rather than the total (set to true if you want primary loop contributors)
     )pbdoc",
         py::arg("tolerance") = 1e-3,
         py::arg("minimumNrOfActiveConnections") = 2,
         py::arg("maximumNrOfActiveConnections") = -1)
    .def("get_nr_of_active_springs",
         &mehp::MEHPForceRelaxation::getNrOfActiveSprings,
         R"pbdoc(
           Get the number of active springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("count_active_clustered_atoms",
         py::overload_cast<const double>(
           &mehp::MEHPForceRelaxation::countActiveClusteredAtoms),
         R"pbdoc(
          Counts the active clustered atoms in the system.

          :param tolerance: springs under this length are considered inactive.)pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_soluble_weight_fraction",
         &mehp::MEHPForceRelaxation::getSolubleWeightFraction,
         R"pbdoc(
          Compute the weight fraction of springs connected to active
          springs (any depth).

          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_dangling_weight_fraction",
         &mehp::MEHPForceRelaxation::getDanglingWeightFraction,
         R"pbdoc(
          Compute the weight fraction of non-active springs

          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_active_chains",
         &mehp::MEHPForceRelaxation::getActiveChains,
         R"pbdoc(
          Get the cross-linker chains that are active.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_effective_functionality_of_atoms",
         &mehp::MEHPForceRelaxation::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
          Returns the number of active springs connected to each atom, atomId used as index

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_spring_lengths",
         &mehp::MEHPForceRelaxation::getSpringLengths,
         R"pbdoc(
          Get the current lengths for all the springs.

          Returns:
               - distances: a vector of size nrOfSprings, with each the norm of the distances
     )pbdoc")
    .def("get_spring_distances",
         &mehp::MEHPForceRelaxation::getSpringDistances,
         R"pbdoc(
          Get the current coordinate differences for all the springs.

          Returns:
               - distances: a vector of size 3*nrOfSprings, with each x, y, z values of the springs
     )pbdoc")
    .def("get_average_spring_length",
         &mehp::MEHPForceRelaxation::getAverageSpringLength,
         R"pbdoc(
           Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()`,
           this value is normalized by the number of springs rather than the number of chains.
     )pbdoc")
    .def("get_default_r0_square",
         &mehp::MEHPForceRelaxation::getDefaultR0Square,
         R"pbdoc(
           Returns the value effectively used in :func:`~pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()` for :math:`\langle R_{0,\eta}^2\rangle`.
     )pbdoc")
    .def("get_nr_of_iterations",
         &mehp::MEHPForceRelaxation::getNrOfIterations,
         R"pbdoc(
          Returns the number of iterations used for force relaxation.
     )pbdoc")
    .def("get_exit_reason",
         &mehp::MEHPForceRelaxation::getExitReason,
         R"pbdoc(
           Returns the reason for termination of the simulation
     )pbdoc")
    .def("requires_another_run",
         &mehp::MEHPForceRelaxation::suggestsRerun,
         R"pbdoc(
          For performance reasons, the objective is only minimised within the distances of one box.
          This means, that there is a possibility, e.g. for a single strand longer than two boxes,
          that it would not be globally minimised.

          If the final displacement of one of the atoms is close
          (1e-3, configurable via :func:`~pylimer_tools_cpp.MEHPForceRelaxation.configRerunEpsilon()`)
          to the imposed min/max, after minimizing,
          this method would return true.
     )pbdoc")
    .def("get_crosslinker_universe",
         &mehp::MEHPForceRelaxation::getCrosslinkerVerse,
         R"pbdoc(
          Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
     )pbdoc")
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const mehp::MEHPForceRelaxation& u) {
        return py::make_tuple(pylimer_tools::utils::serializeToString(u));
      },
      [](py::tuple t) {
        std::string in = t[0].cast<std::string>();
        return mehp::MEHPForceRelaxation::constructFromString(in);
      }))
#endif
    ;

  ////////////////////////////////////////////////////////////////
  // MARK: Configuration Enums
  py::enum_<mehp::StructureSimplificationMode>(m, "StructureSimplificationMode")
    .value("NO_SIMPLIFICATION",
           mehp::StructureSimplificationMode::NO_SIMPLIFICATION)
    .value("X2F_ONLY", mehp::StructureSimplificationMode::X2F_ONLY)
    .value("INACTIVE_ONLY", mehp::StructureSimplificationMode::INACTIVE_ONLY)
    .value("ALL_TIM", mehp::StructureSimplificationMode::ALL_TIM)
    .value("ALL_ANDREI", mehp::StructureSimplificationMode::ALL_ANDREI);

  py::enum_<mehp::LinkSwappingMode>(m, "LinkSwappingMode")
    .value("NO_SWAPPING", mehp::LinkSwappingMode::NO_SWAPPING)
    .value("SLIPLINKS_ONLY", mehp::LinkSwappingMode::SLIPLINKS_ONLY)
    .value("ALL", mehp::LinkSwappingMode::ALL)
    .value("ALL_CYCLE", mehp::LinkSwappingMode::ALL_CYCLE)
    .value("ALL_MC", mehp::LinkSwappingMode::ALL_MC)
    .value("ALL_MC_CYCLE", mehp::LinkSwappingMode::ALL_MC_CYCLE)
    .value("ALL_MC_TRY", mehp::LinkSwappingMode::ALL_MC_TRY)
    .value("ALL_MC_TRY_CYCLE", mehp::LinkSwappingMode::ALL_MC_TRY_CYCLE);

  py::enum_<mehp::SLESolver>(m, "SLESolver")
    .value("DEFAULT", mehp::SLESolver::DEFAULT)
    .value("SIMPLICIAL_LLT", mehp::SLESolver::SIMPLICIAL_LLT)
    .value("SIMPLICIAL_DLT", mehp::SLESolver::SIMPLICIAL_DLT)
    .value("SPARSE_LU", mehp::SLESolver::SPARSE_LU)
    .value("SPARSE_QR", mehp::SLESolver::SPARSE_QR)
    .value("CONJUGATE_GRADIENT", mehp::SLESolver::CONJUGATE_GRADIENT)
    .value("LEAST_SQUARES_CONJUGATE_GRADIENT",
           mehp::SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT)
    .value("BICGSTAB", mehp::SLESolver::BICGSTAB);

  ////////////////////////////////////////////////////////////////
  // MARK: Force Balance
  py::class_<mehp::MEHPForceBalance>(m,
                                     "MEHPForceBalance",
                                     R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.

    This is the second implementation in the group of MEHP provided by this package.
    The distinct feature here is the slip-links: a form of entanglement,
    represented as an entanglement link, just like a four-functional cross-link,
    but with the ability to slip along the two associated strands, therewith
    adjusting the fraction of the contour length on both sides of the link.
     )pbdoc")
    .def(py::init<pe::Universe, int, bool, bool, bool>(),
         R"pbdoc(
          Instantiate the simulator for a certain universe.

          :param universe: the universe to simulate with
          :param crosslinker_type: The atom type of the cross-linkers. Needed to reduce the network.
          :param is2D: Whether to ignore the z direction.
          :param kappa: the spring constant
          :param remove_2functionalCrosslinkers: whether to keep or remove the 2-functional crosslinkers when setting up the network
          :param remove_dangling_chains: whether to keep or remove obviously dangling chains when setting up the network
          )pbdoc",
         py::arg("universe"),
         py::arg("crosslinker_type") = 2,
         py::arg("is_2d") = false,
         py::arg("remove_2functional_crosslinkers") = false,
         py::arg("remove_dangling_chains") = false)
    .def("__copy__",
         [](const mehp::MEHPForceBalance& self) {
           return mehp::MEHPForceBalance(self);
         })
    .def_static("construct_with_random_sliplinks",
                &mehp::MEHPForceBalance::constructWithRandomSlipLinks,
                R"pbdoc(
          Instantiate this simulator with randomly chosen slip-links.
         )pbdoc",
                py::arg("universe"),
                py::arg("nr_of_sliplinks_to_sample"),
                py::arg("upper_sampling_cutoff") = 1.2,
                py::arg("lower_sampling_cutoff") = 0.0,
                py::arg("minimum_nr_of_sliplinks") = 0,
                py::arg("same_strand_cutoff") = 3,
                py::arg("seed") = "",
                py::arg("crosslinker_type") = 2,
                py::arg("is_2d") = false,
                py::arg("skip_dangling_soluble_entanglements") = true)
    .def_property_readonly("network", &mehp::MEHPForceBalance::getNetwork)
    .def(
      "validate_network",
      [](const mehp::MEHPForceBalance& fb) { return fb.validateNetwork(); },
      R"pbdoc(
          Validates the internal structures.

          Throws an error if something is not ok.
          Otherwise, it returns true.

          Can be used e.g. as :code:`assert fb.validate_network()`.
     )pbdoc")
    .def(
      "run_force_relaxation",
      [](mehp::MEHPForceBalance& sim,
         long int maxNrOfSteps,
         // default: 10000
         double xtol,
         const double initialResidualToUse,
         const mehp::StructureSimplificationMode simplificationMode,
         const double inactiveRemovalCutoff,
         bool doInnerIterations,
         const mehp::LinkSwappingMode allowSlipLinksToPassEachOther,
         const int swappingFrequency,
         const double oneOverSpringPartitionUpperLimit,
         const int nrOfCrosslinkSwapsAllowedPerSliplink,
         const bool disableSlipping) {
        return sim.runForceRelaxation(
          maxNrOfSteps,
          xtol,
          initialResidualToUse,
          simplificationMode,
          inactiveRemovalCutoff,
          doInnerIterations,
          allowSlipLinksToPassEachOther,
          swappingFrequency,
          oneOverSpringPartitionUpperLimit,
          nrOfCrosslinkSwapsAllowedPerSliplink,
          disableSlipping,
          []() { return PyErr_CheckSignals() != 0; },
          []() { throw py::error_already_set(); });
      },
      R"pbdoc(
          Run the simulation.
          Note that the final state of the minimization is persisted and reused if you use this method again.
          This is useful if you want to run a global optimization first and add a local one afterwards.
          As a consequence though, you cannot simply benchmark only this method; you must include the setup.

          :param max_nr_of_steps: The maximum number of steps to do during the simulation.
          :param x_tolerance: The tolerance of the displacements as an exit condition.
          :param initial_residual_norm: The residual norm relative to which the relative tolerance is specified. Negative values mean, it will be replaced with the current one.
          :param simplification_mode: How to simplify the structure during the minimization.
          :param inactive_removal_cutoff: The tolerance in distance units of the partial spring length to count as active, relevant if simplification mode is specified to be something other than NO_SIMPLIFICATION.
          :param do_inner_iterations: Whether to do inner iterations; usually, they are not helpful.
          :param allow_sliplinks_to_pass_each_other: Whether slip-links can pass each other.
          :param swapping_frequency: How often slip-links attempt to swap.
          :param one_over_spring_partition_upper_limit: Super-secret parameter. Use 1, gradually increase (and then -1) if you want to publish.
          :param nr_of_crosslink_swaps_allowed_per_sliplink: Use to steer whether slip-links can cross cross-links when swapping is enabled.
          :param disable_slipping: Whether slip-links should be prohibited from slipping.
          )pbdoc",
      py::arg("max_nr_of_steps") = 250000,
      py::arg("x_tolerance") = 1e-12,
      py::arg("initial_residual_norm") = -1.0,
      py::arg("simplification_mode") =
        mehp::StructureSimplificationMode::NO_SIMPLIFICATION,
      py::arg("inactive_removal_cutoff") = 1e-3,
      py::arg("do_inner_iterations") = false,
      py::arg("allow_sliplinks_to_pass_each_other") =
        mehp::LinkSwappingMode::NO_SWAPPING,
      py::arg("swapping_frequency") = 10,
      py::arg("one_over_spring_partition_upper_limit") = 1.0,
      py::arg("nr_of_crosslink_swaps_allowed_per_sliplink") = -1,
      py::arg("disable_slipping") = false)
    .def("deform_to",
         &mehp::MEHPForceBalance::deformTo,
         R"pbdoc(
          Perform a deformation of the system box to a different box.
          All coordinates etc. will be scaled as needed.
         )pbdoc",
         py::arg("new_box"))
    .def("config_step_output",
         &mehp::MEHPForceBalance::configStepOutput,
         R"pbdoc(
          Set which values to log.

          Arguments:
               - values: a list of OutputConfiguration structs
     )pbdoc")
    .def("config_assume_box_large_enough",
         &mehp::MEHPForceBalance::configAssumeBoxLargeEnough,
         R"pbdoc(
          Configure whether to run PBC on the bonds or not.

          If your bonds could get larger than half the box length, this must be kept false (default).
          Otherwise, you can set it to true and therewith get some securities.
         )pbdoc",
         py::arg("box_large_enough") = false)
    .def("config_mean_bond_length",
         &mehp::MEHPForceBalance::configMeanBondLength,
         R"pbdoc(
     Configure the :math:`b` used e.g. for the topological Gamma-factor.
     )pbdoc",
         py::arg("b") = 1.0)
    .def("config_spring_constant",
         &mehp::MEHPForceBalance::configSpringConstant,
         R"pbdoc()pbdoc",
         py::arg("kappa") = 1.0)
    .def("config_entanglement_type",
         &mehp::MEHPForceBalance::configEntanglementType,
         R"pbdoc(
         To have certain cross-links behave as entanglements in the removal process,
         you can specify the here a type, that you have used in the universe to specify:
         - the type of entanglement atoms (expected with functionality f = 3),
         - and the entanglement-bonds between the entanglement atoms.

         I.e., say you want to model some entanglements as non-slipping,
         bonds between two strand beads resulting in f = 3 beads, for example,
         you can call this method to have the "StructureSimplificationMode" also remove these atoms,
         if they have a functionality of 2 or less while still being connected to its partner bead.
         )pbdoc",
         py::arg("type") = -1)
    .def("config_spring_breaking_distance",
         &mehp::MEHPForceBalance::configSpringBreakingDistance,
         R"pbdoc(
          Configure the "force" (distance over contour length) at which the bonds break.
          Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
          Springs breaking will happen before the simplification procedure is run.
          Negative values will disable spring breaking.
          Default: -1..
         )pbdoc",
         py::arg("distance_over_contour_length") = -1)
    .def("config_simplification_frequency",
         &mehp::MEHPForceBalance::configSimplificationFrequency,
         R"pbdoc(
         Config every how many steps to simplify the structure.
         Default: 10.
         )pbdoc",
         py::arg("frequency") = 10)
    .def("swap_sliplinks_incl_xlinks",
         &mehp::MEHPForceBalance::swapSlipLinksInclXlinks)
    .def("move_sliplinks_to_their_best_branch",
         &mehp::MEHPForceBalance::moveSlipLinksToTheirBestBranch)
    .def(
      "get_force_on",
      [](mehp::MEHPForceBalance& sim,
         const size_t linkIdx,
         const double oneOver) { return sim.getForceOn(linkIdx, oneOver); },
      R"pbdoc(
          Evaluate the force on a particular (slip- or cross-) link.
      )pbdoc",
      py::arg("link_idx"),
      py::arg("one_over_spring_partition_upper_limit") = 1.0)
    .def("get_force_magnitude_vector",
         &mehp::MEHPForceBalance::getForceMagnitudeVector,
         R"pbdoc(
          Evaluate the norm of the force on each (slip- or cross-) link.
     )pbdoc")
    .def(
      "get_stress_on",
      [](mehp::MEHPForceBalance& sim,
         const size_t linkIdx,
         const double oneOver) { return sim.getStressOn(linkIdx, oneOver); },
      R"pbdoc(
          Evaluate the stress on a particular (slip- or cross-) link.
      )pbdoc",
      py::arg("link_idx"),
      py::arg("one_over_spring_partition_upper_limit") = 1.0)
    .def("inspect_displacement_to_mean_position_update",
         &mehp::MEHPForceBalance::inspectDisplacementToMeanPositionUpdate,
         R"pbdoc(
          Helper method to debug and/or understand what happens to certain links when being displaced.
         )pbdoc",
         py::arg("link_idx"),
         py::arg("one_over_spring_partition_upper_limit") = 1.0)
    .def("inspect_spring_partition_update",
         &mehp::MEHPForceBalance::inspectSpringPartitionUpdate,
         R"pbdoc(
          Helper method to debug and/or understand what happens to certain links
          when the spring partition is being updated.
         )pbdoc",
         py::arg("link_idx"))
    .def("inspect_parametrisation_optimsation_for_link",
         &mehp::MEHPForceBalance::inspectParametrisationOptimsationForLink,
         R"pbdoc(
          Helper method to debug and/or understand what happens to certain links
          when being displaced and their partition updated.
         )pbdoc",
         py::arg("link_idx"),
         py::arg("displacements"),
         py::arg("spring_partitions"),
         py::arg("max_nr_of_steps") = 100,
         py::arg("alpha_tol") = 1e-9,
         py::arg("min_nr_of_steps") = 1,
         py::arg("one_over_spring_partition_upper_limit") = 1.0)
    .def("get_springpartition_indices_of_sliplink",
         &mehp::MEHPForceBalance::getSpringpartitionIndicesOfSliplink,
         R"pbdoc()pbdoc",
         py::arg("network"),
         py::arg("link_idx"))
    .def("get_neighbour_link_indices",
         &mehp::MEHPForceBalance::getNeighbourLinkIndices,
         R"pbdoc()pbdoc",
         py::arg("network"),
         py::arg("link_idx"))
    .def(
      "evaluate_partial_spring_distance",
      [](const mehp::MEHPForceBalance& sim,
         const mehp::ForceBalanceNetwork& net,
         const Eigen::VectorXd& u,
         const size_t springIdx) {
        return sim.evaluatePartialSpringDistance(net, u, springIdx);
      },
      R"pbdoc()pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"))
    .def(
      "evaluate_partial_spring_distance_from",
      [](const mehp::MEHPForceBalance& sim,
         const mehp::ForceBalanceNetwork& net,
         const Eigen::VectorXd& u,
         const size_t springIdx,
         const size_t linkIdx) {
        return sim.evaluatePartialSpringDistanceFrom(
          net, u, springIdx, linkIdx);
      },
      R"pbdoc()pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"),
      py::arg("link_idx"))
    .def(
      "evaluate_partial_spring_distance_to",
      [](const mehp::MEHPForceBalance& sim,
         const mehp::ForceBalanceNetwork& net,
         const Eigen::VectorXd& u,
         const size_t springIdx,
         const size_t linkIdx) {
        return sim.evaluatePartialSpringDistanceTo(net, u, springIdx, linkIdx);
      },
      R"pbdoc()pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"),
      py::arg("link_idx"))
    // .def("getForceEvaluator", &mehp::MEHPForceBalance::getForceEvaluator,
    // R"pbdoc(
    //      Query the currently used force evaluator.
    // )pbdoc")
    //     .def("setForceEvaluator",
    //          &mehp::MEHPForceBalance::setForceEvaluator,
    //          R"pbdoc(
    //           Reset the currently used force evaluator.
    //      )pbdoc")
    //     .def("getForce",
    //          &mehp::MEHPForceBalance::getForce,
    //          R"pbdoc(
    //           Returns the force at the current state of the simulation.
    //      )pbdoc")
    //     .def("getResidualNorm",
    //          &mehp::MEHPForceBalance::getResidualNorm,
    //          R"pbdoc(
    //           Returns the residual norm at the current state of the
    //           simulation.
    //      )pbdoc")
    .def("get_nr_of_atoms", &mehp::MEHPForceBalance::getNumAtoms)
    .def("get_nr_of_extra_atoms", &mehp::MEHPForceBalance::getNumExtraAtoms)
    .def("get_nr_of_bonds", &mehp::MEHPForceBalance::getNumBonds)
    .def("get_nr_of_extra_bonds", &mehp::MEHPForceBalance::getNumExtraBonds)
    .def("get_nr_of_intra_chain_sliplinks",
         &mehp::MEHPForceBalance::getNumIntraChainSlipLinks)
    .def("get_pressure",
         &mehp::MEHPForceBalance::getPressure,
         R"pbdoc(
          Returns the pressure at the current state of the simulation.
     )pbdoc")
    .def("get_soluble_weight_fraction",
         &mehp::MEHPForceBalance::getSolubleWeightFraction,
         R"pbdoc(
          Compute the weight fraction of springs connected to active
          springs (any depth).

          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_dangling_weight_fraction",
         &mehp::MEHPForceBalance::getDanglingWeightFraction,
         R"pbdoc(
          Compute the weight fraction of non-active springs

          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("add_sliplinks",
         py::overload_cast<const std::vector<size_t>&,
                           const std::vector<size_t>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const bool>(&mehp::MEHPForceBalance::addSlipLinks),
         R"pbdoc(
          Add new slip-links
     )pbdoc",
         py::arg("strand_idx_1"),
         py::arg("strand_idx_2"),
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("alpha_1"),
         py::arg("alpha_2"),
         py::arg("clamp_alpha") = false)
    .def("randomly_add_sliplinks",
         &mehp::MEHPForceBalance::randomlyAddSliplinks,
         R"pbdoc(
          Randomly sample and add slip-links based on certain criteria.
         )pbdoc",
         py::arg("nr_of_sliplinks_to_sample"),
         py::arg("cutoff") = 2.0,
         py::arg("minimum_nr_of_sliplinks") = 0,
         py::arg("same_strand_cutoff") = 2,
         py::arg("exclude_crosslinks") = false,
         py::arg("seed") = -1)
    .def("add_sliplinks_based_on_cycles",
         &mehp::MEHPForceBalance::addSliplinksBasedOnCycles,
         R"pbdoc(
          Detect and add slip-links based on detected entanglements.

          WARNING:
               Does not work yet.
         )pbdoc",
         py::arg("maxLoopLength") = -1)
    .def(
      "get_stress_tensor",
      [](mehp::MEHPForceBalance& fb, const double oneOver = 1.) {
        return fb.getStressTensor(oneOver);
      },
      R"pbdoc(
          Returns the stress tensor at the current state of the simulation.

          The units are in :math:`[\text{units of }\kappa]/[\text{distance units}]`,
          where the units of :math:`\kappa` should be :math:`[\text{force}]/[\text{distance units}]^2`.
          Make sure to multiply by :math:`\kappa` or configure it appropriately.
     )pbdoc",
      py::arg("one_over_spring_partition_upper_limit") = 1.)
    .def("get_stress_tensor_link_based",
         &mehp::MEHPForceBalance::getStressTensorLinkBased,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc",
         py::arg("one_over_spring_partition_upper_limit") = 1.,
         py::arg("xlinks_only") = false)
    .def("get_gamma_factor",
         &mehp::MEHPForceBalance::getGammaFactor,
         R"pbdoc(
          Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

          :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
          which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
          where :math:`\eta` is the index of a particular strand,
          :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\eta}*b^2$`
          :math:`N_{\eta}` is the number of atoms in this strand :math:`\eta`,
          :math:`b` its mean square bond length,
          :math:`T` the temperature and
          :math:`k_B` Boltzmann's constant.

          :param b02: the melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\infinity b^2`.
          :param nr_of_chains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of springs.
     )pbdoc",
         py::arg("b02") = -1.0,
         py::arg("nr_of_chains") = -1,
         py::arg("one_over_spring_partition_upper_limit") = 1.)
    .def("get_gamma_factors",
         &mehp::MEHPForceBalance::getGammaFactors,
         R"pbdoc(
          Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)
     )pbdoc",
         py::arg("b02"),
         py::arg("one_over_spring_partition_upper_limit") = 1.)
    .def("get_gamma_factors_in_dir",
         &mehp::MEHPForceBalance::getGammaFactorsInDir,
         R"pbdoc(
               Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)

               :param b02: the melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\infinity b^2`.
               :param direction: the direction in which to compute the gamma factors (0: x, 1: y, 2: z)
          )pbdoc",
         py::arg("b02"),
         py::arg("direction"),
         py::arg("one_over_spring_partition_upper_limit") = 1.)
    .def("get_nr_of_nodes",
         &mehp::MEHPForceBalance::getNrOfNodes,
         R"pbdoc(
           Get the number of nodes (crosslinkers) considered in this simulation.
     )pbdoc")
    .def("get_nr_of_springs",
         &mehp::MEHPForceBalance::getNrOfSprings,
         R"pbdoc(
          Get the number of springs considered in this simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc")
    .def("get_spring_partitions",
         &mehp::MEHPForceBalance::getSpringPartitions,
         R"pbdoc(
          Get the current spring partitions (the fraction of the contour length associated with each partial spring).
     )pbdoc")
    .def("get_weighted_partial_spring_lengths",
         &mehp::MEHPForceBalance::getWeightedPartialSpringLengths,
         R"pbdoc(
          Get the current partial spring lengths (norm of vector) divided by the spring partition times the contour length.
          )pbdoc",
         py::arg("one_over_spring_partition_upper_limit") = 1.)
    .def("set_spring_partitions",
         &mehp::MEHPForceBalance::setSpringPartitions,
         R"pbdoc(
          Set the current spring partitions.
     )pbdoc")
    .def("get_displacements",
         &mehp::MEHPForceBalance::getCurrentDisplacements,
         R"pbdoc(
          Get the current link displacements.
     )pbdoc")
    .def("set_displacements",
         &mehp::MEHPForceBalance::setCurrentDisplacements,
         R"pbdoc(
          Set the current link displacements.
     )pbdoc")
    .def("set_spring_contour_lengths",
         &mehp::MEHPForceBalance::setSpringContourLengths,
         R"pbdoc(
          Set/overwrite the contour lengths.
     )pbdoc")
    .def("get_displacement_residual_norm",
         &mehp::MEHPForceBalance::getDisplacementResidualNorm,
         R"pbdoc(
          Get the current link displacement residual norm.
     )pbdoc",
         py::arg("one_over_spring_partition_upper_limit") = 1.)
    .def("get_ids_of_active_nodes",
         &mehp::MEHPForceBalance::getIdsOfActiveNodes,
         R"pbdoc(
          Get the atom ids of the nodes that are considered active.
          Only cross-link ids are returned (not e.g. entanglement links).

          :param tolerance: springs under this length are considered inactive. A node is active if it has > 1 active springs.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_nodes",
         &mehp::MEHPForceBalance::getNrOfActiveNodes,
         R"pbdoc(
          Get the number of active nodes (incl. entanglement nodes [atoms with type = entanglementType, present in the universe when creating this simulator],
          excl. entanglement links [the slip-links created internally when e.g. constructing the simulator with random slip-links]).

          :param tolerance: springs under this length are considered inactive. A node is active if it has > 1 active springs.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_springs",
         &mehp::MEHPForceBalance::getNrOfActiveSprings,
         R"pbdoc(
           Get the number of active springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_springs_in_dir",
         &mehp::MEHPForceBalance::getNrOfActiveSpringsInDir,
         R"pbdoc(
                Get the number of active springs remaining after running the simulation.

               :param direction: the direction in which to compute the active springs (0: x, 1: y, 2: z)
               :param tolerance: springs under this length are considered inactive
          )pbdoc",
         py::arg("direction"),
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_partial_springs",
         &mehp::MEHPForceBalance::getNrOfActivePartialSprings,
         R"pbdoc(
           Get the number of active partial springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_current_partial_spring_vectors",
         &mehp::MEHPForceBalance::getCurrentPartialSpringDistances,
         R"pbdoc(
          Get the partial spring vectors.
     )pbdoc")
    .def("get_current_partial_spring_lengths",
         &mehp::MEHPForceBalance::getCurrentPartialSpringLengths,
         R"pbdoc(
          Get the partial spring distances.
     )pbdoc")
    .def("get_current_spring_vectors",
         &mehp::MEHPForceBalance::getCurrentSpringDistances)
    .def("get_overall_spring_lengths",
         &mehp::MEHPForceBalance::getOverallSpringLengths,
         R"pbdoc(
          Get the sum of the lengths of the partial springs of each spring.
     )pbdoc")
    .def("get_effective_functionality_of_atoms",
         &mehp::MEHPForceBalance::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
          Returns the number of active springs connected to each atom, atomId used as index

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_average_spring_length",
         &mehp::MEHPForceBalance::getAverageSpringLength,
         R"pbdoc(
           Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()`,
           this value is normalized by the number of springs rather than the number of chains.
     )pbdoc")
    .def("get_default_mean_bond_length",
         &mehp::MEHPForceBalance::getDefaultMeanBondLength,
         R"pbdoc(
           Returns the value effectively used in :func:`~pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()` for
           :math:`b` in :math:`\langle R_{0,\eta}^2 = N_{\eta} b^2\rangle`.
     )pbdoc")
    .def("get_nr_of_iterations",
         &mehp::MEHPForceBalance::getNrOfIterations,
         R"pbdoc(
          Returns the number of iterations used for force relaxation so far.
     )pbdoc")
    .def("get_exit_reason",
         &mehp::MEHPForceBalance::getExitReason,
         R"pbdoc(
           Returns the reason for termination of the simulation
     )pbdoc")
    .def("get_crosslinker_universe",
         &mehp::MEHPForceBalance::getCrosslinkerVerse,
         R"pbdoc(
          Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
     )pbdoc")
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const mehp::MEHPForceBalance& u) {
        return py::make_tuple(pylimer_tools::utils::serializeToString(u));
      },
      [](py::tuple t) {
        std::string in = t[0].cast<std::string>();
        return mehp::MEHPForceBalance::constructFromString(in);
      }))
#endif
    ;

  ////////////////////////////////////////////////////////////////
  // MARK: Force Balance 2
  py::class_<mehp::MEHPForceBalance2>(m,
                                      "MEHPForceBalance2",
                                      R"pbdoc(
     A small simulation tool for quickly minimizing the force between the cross-linker beads.

     This is the third implementation of the MEHP. 
     It's the fastest implementation by using a simple spring model only, disabling the non-linear
     periodic boundary conditions, and instead builds a sparse linear system of equations that's readily solved.
     However, it allows entanglements to be represented as additional links or/and springs,
     although without slipping along the chain.
      )pbdoc")
    .def(py::init<pe::Universe, int, bool>(),
         R"pbdoc(
           Instantiate the simulator for a certain universe.

           :param universe: the universe giving the basic connectivity to compute with
           :param crosslinker_type: The atom type of the cross-linkers. Needed to reduce the network.
           :param is_2d: Whether to ignore the z direction.
           :param kappa: the spring constant
           :param remove_2functionalCrosslinkers: whether to keep or remove the 2-functional crosslinkers when setting up the network
           :param remove_dangling_chains: whether to keep or remove obviously dangling chains when setting up the network
           )pbdoc",
         py::arg("universe"),
         py::arg("crosslinker_type") = 2,
         py::arg("is_2d") = false)
    .def(py::init<
           pe::Universe,
           pylimer_tools::topo::entanglement_detection::AtomPairEntanglements,
           int,
           bool,
           bool>(),
         R"pbdoc(
           Instantiate the simulator for a certain universe with the given entanglements.

           :param universe: the universe giving the basic connectivity to compute with
           :param entanglements: the entanglements to use in the computation
           :param crosslinker_type: The atom type of the cross-linkers. Needed to reduce the network.
           :param is_2d: Whether to ignore the z direction.
           :param entanglements_as_springs: whether to use the entanglements as springs instead of links
           )pbdoc",
         py::arg("universe"),
         py::arg("entanglements"),
         py::arg("crosslinker_type") = 2,
         py::arg("is_2d") = false,
         py::arg("entanglements_as_springs") = false)

    .def(py::init<pe::Universe,
                  size_t,
                  double,
                  double,
                  size_t,
                  double,
                  std::string,
                  int,
                  bool,
                  bool,
                  bool>(),
         R"pbdoc(
     Instantiate this simulator with randomly chosen slip-links.

     :param universe: the universe containing the basic atoms and connectivity
     :param nr_of_entanglements_to_sample: the number of entanglements to sample
     :param upper_cutoff: maximum distance from one sampled bead to its partner
     :param lower_cutoff: minimum distance from one sampled bead to its partner
     :param minimum_nr_of_entanglements: the minimum number of entanglements that should be sampled
     :param same_strand_cutoff: distance from one sampled bead to its pair within the same strand
     :param seed: the seed for the random number generator
     :param cross_linker_type:
     :param is_2d:
     :param filter_entanglements:
     :param entanglements_as_springs: whether to model the entanglements as merged beads or beads with 1 spring in between
   )pbdoc",
         py::arg("universe"),
         py::arg("nr_of_entanglements_to_sample"),
         py::arg("upper_sampling_cutoff") = 1.2,
         py::arg("lower_sampling_cutoff") = 0.0,
         py::arg("minimum_nr_of_sliplinks") = 0,
         py::arg("same_strand_cutoff") = 3,
         py::arg("seed") = "",
         py::arg("crosslinker_type") = 2,
         py::arg("is_2d") = false,
         py::arg("skip_dangling_soluble_entanglements") = true,
         py::arg("entanglements_as_springs") = true)

    .def("__copy__",
         [](const mehp::MEHPForceBalance2& self) {
           return mehp::MEHPForceBalance2(self);
         })
    .def_property_readonly("network", &mehp::MEHPForceBalance2::getNetwork)
    .def(
      "validate_network",
      [](const mehp::MEHPForceBalance2& fb) { return fb.validateNetwork(); },
      R"pbdoc(
           Validates the internal structures.

           Throws an error if something is not ok.
           Otherwise, it returns true.

           Can be used e.g. as :code:`assert fb.validate_network()`.
      )pbdoc")
    .def(
      "run_force_relaxation",
      [](mehp::MEHPForceBalance2& sim,
         const mehp::StructureSimplificationMode simplificationMode,
         const double inactiveRemovalCutoff,
         const mehp::SLESolver solverChoice) {
        return sim.runForceRelaxation(
          simplificationMode,
          inactiveRemovalCutoff,
          solverChoice,
          []() { return PyErr_CheckSignals() != 0; },
          []() { throw py::error_already_set(); });
      },
      R"pbdoc(
           Run the simulation.

           :param simplification_mode: How to simplify the structure during the minimization.
           :param inactive_removal_cutoff: The tolerance in distance units of the partial spring length to count as active, relevant if simplification mode is specified to be something other than NO_SIMPLIFICATION.
           :param sle_solver: The solver to use for the system of linear equations (SLE).
           )pbdoc",
      py::arg("simplification_mode") =
        mehp::StructureSimplificationMode::NO_SIMPLIFICATION,
      py::arg("inactive_removal_cutoff") = 1e-3,
      py::arg("sle_solver") = mehp::SLESolver::DEFAULT)
    .def("deform_to",
         &mehp::MEHPForceBalance2::deformTo,
         R"pbdoc(
           Perform a deformation of the system box to a different box.
           All coordinates etc. will be scaled as needed.
          )pbdoc",
         py::arg("new_box"))
    .def("config_step_output",
         &mehp::MEHPForceBalance2::configStepOutput,
         R"pbdoc(
           Set which values to log.

           Arguments:
                - values: a list of OutputConfiguration structs
      )pbdoc")
    .def("config_mean_bond_length",
         &mehp::MEHPForceBalance2::configMeanBondLength,
         R"pbdoc(
      Configure the :math:`b` used e.g. for the topological Gamma-factor.
      )pbdoc",
         py::arg("b") = 1.0)
    .def("config_spring_constant",
         &mehp::MEHPForceBalance2::configSpringConstant,
         R"pbdoc()pbdoc",
         py::arg("kappa") = 1.0)
    .def("config_spring_breaking_distance",
         &mehp::MEHPForceBalance2::configSpringBreakingDistance,
         R"pbdoc(
           Configure the "force" (distance over contour length) at which the bonds break.
           Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
           Springs breaking will happen before the simplification procedure is run.
           Negative values will disable spring breaking.
           Default: -1..
          )pbdoc",
         py::arg("distance_over_contour_length") = -1)

    .def(
      "get_force_on",
      [](mehp::MEHPForceBalance2& sim, const size_t linkIdx) {
        return sim.getForceOn(linkIdx);
      },
      R"pbdoc(
           Evaluate the force on a particular (slip- or cross-) link.
       )pbdoc",
      py::arg("link_idx"))
    .def("get_force_magnitude_vector",
         &mehp::MEHPForceBalance2::getForceMagnitudeVector,
         R"pbdoc(
           Evaluate the norm of the force on each (slip- or cross-) link.
      )pbdoc")
    .def(
      "get_stress_on",
      [](mehp::MEHPForceBalance2& sim, const size_t linkIdx) {
        return sim.getStressOn(linkIdx);
      },
      R"pbdoc(
           Evaluate the stress on a particular (slip- or cross-) link.
       )pbdoc",
      py::arg("link_idx"))
    .def("inspect_displacement_to_mean_position_update",
         &mehp::MEHPForceBalance2::inspectDisplacementToMeanPositionUpdate,
         R"pbdoc(
           Helper method to debug and/or understand what happens to certain links when being displaced.
          )pbdoc",
         py::arg("link_idx"))
    .def("get_neighbour_link_indices",
         &mehp::MEHPForceBalance2::getNeighbourLinkIndices,
         R"pbdoc()pbdoc",
         py::arg("network"),
         py::arg("link_idx"))
    .def(
      "evaluate_partial_spring_distance",
      [](const mehp::MEHPForceBalance2& sim,
         const mehp::ForceBalance2Network& net,
         const Eigen::VectorXd& u,
         const size_t springIdx) {
        return sim.evaluateSpringVector(net, u, springIdx);
      },
      R"pbdoc()pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"))
    .def(
      "evaluate_partial_spring_distance_from",
      [](const mehp::MEHPForceBalance2& sim,
         const mehp::ForceBalance2Network& net,
         const Eigen::VectorXd& u,
         const size_t springIdx,
         const size_t linkIdx) {
        return sim.evaluateSpringVectorFrom(net, u, springIdx, linkIdx);
      },
      R"pbdoc()pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"),
      py::arg("link_idx"))
    .def(
      "evaluate_partial_spring_distance_to",
      [](const mehp::MEHPForceBalance2& sim,
         const mehp::ForceBalance2Network& net,
         const Eigen::VectorXd& u,
         const size_t springIdx,
         const size_t linkIdx) {
        return sim.evaluateSpringVectorTo(net, u, springIdx, linkIdx);
      },
      R"pbdoc()pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"),
      py::arg("link_idx"))
    // .def("getForceEvaluator", &mehp::MEHPForceBalance2::getForceEvaluator,
    // R"pbdoc(
    //      Query the currently used force evaluator.
    // )pbdoc")
    //     .def("setForceEvaluator",
    //          &mehp::MEHPForceBalance2::setForceEvaluator,
    //          R"pbdoc(
    //           Reset the currently used force evaluator.
    //      )pbdoc")
    //     .def("getForce",
    //          &mehp::MEHPForceBalance2::getForce,
    //          R"pbdoc(
    //           Returns the force at the current state of the simulation.
    //      )pbdoc")
    //     .def("getResidualNorm",
    //          &mehp::MEHPForceBalance2::getResidualNorm,
    //          R"pbdoc(
    //           Returns the residual norm at the current state of the
    //           simulation.
    //      )pbdoc")
    .def("get_nr_of_atoms", &mehp::MEHPForceBalance2::getNumAtoms)
    .def("get_nr_of_extra_atoms", &mehp::MEHPForceBalance2::getNumExtraAtoms)
    .def("get_nr_of_bonds", &mehp::MEHPForceBalance2::getNumBonds)
    .def("get_nr_of_extra_bonds", &mehp::MEHPForceBalance2::getNumExtraBonds)
    .def("get_nr_of_intra_chain_sliplinks",
         &mehp::MEHPForceBalance2::getNumIntraChainSlipLinks)
    .def("get_pressure",
         &mehp::MEHPForceBalance2::getPressure,
         R"pbdoc(
           Returns the pressure at the current state of the simulation.
      )pbdoc")
    .def("get_soluble_weight_fraction",
         &mehp::MEHPForceBalance2::getSolubleWeightFraction,
         R"pbdoc(
           Compute the weight fraction of springs connected to active
           springs (any depth).

           Caution: ignores atom masses.
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_dangling_weight_fraction",
         &mehp::MEHPForceBalance2::getDanglingWeightFraction,
         R"pbdoc(
           Compute the weight fraction of non-active springs

           Caution: ignores atom masses.
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def(
      "get_stress_tensor",
      [](mehp::MEHPForceBalance2& fb) { return fb.getStressTensor(); },
      R"pbdoc(
           Returns the stress tensor at the current state of the simulation.

           The units are in :math:`[\text{units of }\kappa]/[\text{distance units}]`,
           where the units of :math:`\kappa` should be :math:`[\text{force}]/[\text{distance units}]^2`.
           Make sure to multiply by :math:`\kappa` or configure it appropriately.
      )pbdoc")
    .def("get_stress_tensor_link_based",
         &mehp::MEHPForceBalance2::getStressTensorLinkBased,
         R"pbdoc(
           Returns the stress tensor at the current state of the simulation.
      )pbdoc",
         py::arg("xlinks_only") = false)
    .def("get_gamma_factor",
         &mehp::MEHPForceBalance2::getGammaFactor,
         R"pbdoc(
           Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

           :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
           which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
           where :math:`\eta` is the index of a particular strand,
           :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\eta}*b^2$`
           :math:`N_{\eta}` is the number of atoms in this strand :math:`\eta`,
           :math:`b` its mean square bond length,
           :math:`T` the temperature and
           :math:`k_B` Boltzmann's constant.

           :param b02: the melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\infinity b^2`.
           :param nr_of_chains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of springs.
      )pbdoc",
         py::arg("b02") = -1.0,
         py::arg("nr_of_chains") = -1)
    .def("get_gamma_factors",
         &mehp::MEHPForceBalance2::getGammaFactors,
         R"pbdoc(
           Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)
      )pbdoc",
         py::arg("b02"))
    .def("get_gamma_factors_in_dir",
         &mehp::MEHPForceBalance2::getGammaFactorsInDir,
         R"pbdoc(
                Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)

                :param b02: the melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\infinity b^2`.
                :param direction: the direction in which to compute the gamma factors (0: x, 1: y, 2: z)
           )pbdoc",
         py::arg("b02"),
         py::arg("direction"))
    .def("get_nr_of_nodes",
         &mehp::MEHPForceBalance2::getNrOfNodes,
         R"pbdoc(
            Get the number of nodes (crosslinkers) considered in this simulation.
      )pbdoc")
    .def("get_nr_of_springs",
         &mehp::MEHPForceBalance2::getNrOfStrands,
         R"pbdoc(
           Get the number of springs considered in this simulation.

           :param tolerance: springs under this length are considered inactive
      )pbdoc")
    .def("get_weighted_partial_spring_lengths",
         &mehp::MEHPForceBalance2::getWeightedSpringLengths,
         R"pbdoc(
           Get the current partial spring lengths (norm of vector) divided by the spring partition times the contour length.
           )pbdoc")
    .def("get_displacements",
         &mehp::MEHPForceBalance2::getCurrentDisplacements,
         R"pbdoc(
           Get the current link displacements.
      )pbdoc")
    .def("set_displacements",
         &mehp::MEHPForceBalance2::setCurrentDisplacements,
         R"pbdoc(
           Set the current link displacements.
      )pbdoc")
    .def("set_spring_contour_lengths",
         &mehp::MEHPForceBalance2::setSpringContourLengths,
         R"pbdoc(
           Set/overwrite the contour lengths.
      )pbdoc")
    .def("get_displacement_residual_norm",
         &mehp::MEHPForceBalance2::getDisplacementResidualNorm,
         R"pbdoc(
           Get the current link displacement residual norm.
      )pbdoc")
    .def("get_ids_of_active_nodes",
         &mehp::MEHPForceBalance2::getIdsOfActiveNodes,
         R"pbdoc(
           Get the atom ids of the nodes that are considered active.
           Only cross-link ids are returned (not e.g. entanglement links).

           :param tolerance: springs under this length are considered inactive. A node is active if it has > 1 active springs.
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_nodes",
         &mehp::MEHPForceBalance2::getNrOfActiveNodes,
         R"pbdoc(
           Get the number of active nodes (incl. entanglement nodes [atoms with type = entanglementType, present in the universe when creating this simulator],
           excl. entanglement links [the slip-links created internally when e.g. constructing the simulator with random slip-links]).

           :param tolerance: springs under this length are considered inactive. A node is active if it has > 1 active springs.
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_springs",
         &mehp::MEHPForceBalance2::getNrOfActiveStrands,
         R"pbdoc(
            Get the number of active springs remaining after running the simulation.

           :param tolerance: springs under this length are considered inactive
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_springs_in_dir",
         &mehp::MEHPForceBalance2::getNrOfActiveStrandsInDir,
         R"pbdoc(
                 Get the number of active springs remaining after running the simulation.

                :param direction: the direction in which to compute the active springs (0: x, 1: y, 2: z)
                :param tolerance: springs under this length are considered inactive
           )pbdoc",
         py::arg("direction"),
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_partial_springs",
         &mehp::MEHPForceBalance2::getNrOfActiveSprings,
         R"pbdoc(
            Get the number of active partial springs remaining after running the simulation.

           :param tolerance: springs under this length are considered inactive
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_current_partial_spring_vectors",
         &mehp::MEHPForceBalance2::getCurrentSpringDistances,
         R"pbdoc(
           Get the partial spring vectors.
      )pbdoc")
    .def("get_current_partial_spring_lengths",
         &mehp::MEHPForceBalance2::getCurrentSpringLengths,
         R"pbdoc(
           Get the partial spring distances.
      )pbdoc")
    .def("get_overall_spring_lengths",
         &mehp::MEHPForceBalance2::getOverallSpringLengths,
         R"pbdoc(
           Get the sum of the lengths of the partial springs of each spring.
      )pbdoc")
    .def("get_effective_functionality_of_atoms",
         &mehp::MEHPForceBalance2::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
           Returns the number of active springs connected to each atom, atomId used as index

           :param tolerance: springs under this length are considered inactive
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_average_spring_length",
         &mehp::MEHPForceBalance2::getAverageSpringLength,
         R"pbdoc(
            Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()`,
            this value is normalized by the number of springs rather than the number of chains.
      )pbdoc")
    .def("get_default_mean_bond_length",
         &mehp::MEHPForceBalance2::getDefaultMeanBondLength,
         R"pbdoc(
            Returns the value effectively used in :func:`~pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()` for
            :math:`b` in :math:`\langle R_{0,\eta}^2 = N_{\eta} b^2\rangle`.
      )pbdoc")
    .def("get_nr_of_iterations",
         &mehp::MEHPForceBalance2::getNrOfIterations,
         R"pbdoc(
           Returns the number of iterations used for force relaxation so far.
      )pbdoc")
    .def("get_exit_reason",
         &mehp::MEHPForceBalance2::getExitReason,
         R"pbdoc(
            Returns the reason for termination of the simulation
      )pbdoc")
    .def("get_crosslinker_universe",
         &mehp::MEHPForceBalance2::getCrosslinkerVerse,
         R"pbdoc(
           Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
      )pbdoc")
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const mehp::MEHPForceBalance2& u) {
        return py::make_tuple(pylimer_tools::utils::serializeToString(u));
      },
      [](py::tuple t) {
        std::string in = t[0].cast<std::string>();
        return mehp::MEHPForceBalance2::constructFromString(in);
      }))
#endif
    ;
  /**
   * DPD Simulations
   */

  ////////////////////////////////////////////////////////////////
  // MARK: DPD Simulator
  py::class_<dpd::DPDSimulator>(m,
                                "DPDSimulator",
                                R"pbdoc(
          A quick-and-dirty implementation of the DPD simulation
          with slip-springs as presented by Langeloth et al.
     )pbdoc")
    .def(py::init<const pe::Universe,
                  const int,
                  const int,
                  const bool,
                  const std::string>(),
         "Get an instance of this class",
         py::arg("universe"),
         py::arg("crosslinker_type") = 2,
         py::arg("slipspring_bond_type") = 9,
         py::arg("is_2d") = false,
         py::arg("seed") = "")
    //     .def("runSimulation",
    //          &dpd::DPDSimulator::runSimulation,
    //          R"pbdoc(
    //           Actually do some simulation steps.
    //      )pbdoc",
    //          py::arg("n_steps"),
    //          py::arg("dt") = 0.06,
    //          py::arg("with_MC") = false)
    .def(
      "run_simulation",
      [](dpd::DPDSimulator& sim, int nSteps, double dt, bool withMC) {
        sim.configTimeStep(dt);
        return sim.runSimulation(
          nSteps,
          withMC,
          []() { return PyErr_CheckSignals() != 0; },
          []() { throw py::error_already_set(); });
      },
      py::arg("n_steps"),
      py::arg("dt") = 0.06,
      py::arg("with_MC") = false)
    .def("assume_box_large_enough",
         &dpd::DPDSimulator::configAssumeBoxLargeEnough,
         R"pbdoc(
          Configure whether to run PBC on the bonds or not.

          If your bonds could get larger than half the box length, this must be kept false (default).
          Otherwise, you can set it to true and therewith get some securities.
         )pbdoc")
    .def("create_slip_springs",
         &dpd::DPDSimulator::createSlipSprings,
         R"pbdoc(
          Randomly add the specified number of slip-springs to neighbours within the specified cut-offs.
     )pbdoc",
         py::arg("num"),
         py::arg("bond_type") = 9)
    .def("config_a",
         &dpd::DPDSimulator::configA,
         R"pbdoc(
          Configure the force-field (pair-style) parameter `A`.
     )pbdoc",
         py::arg("A") = 25.)
    .def("config_sigma",
         &dpd::DPDSimulator::configSigma,
         R"pbdoc(
          Configure the force-field (pair-style) parameter `\sigma`.
     )pbdoc",
         py::arg("sigma") = 3.)
    .def("config_spring_constant",
         &dpd::DPDSimulator::configSpringConstant,
         R"pbdoc(
          Configure the force-field (bond-style) parameter `k`, the spring constant.
     )pbdoc",
         py::arg("k") = 2.)
    .def("config_lambda",
         &dpd::DPDSimulator::configLambda,
         R"pbdoc(
          Configure the modified velocity verlet integration parameter `\lambda`.
     )pbdoc",
         py::arg("l") = 0.65)
    .def("config_slipspring_high_cutoff",
         &dpd::DPDSimulator::configSlipspringHighCutoff,
         R"pbdoc(
          Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.
     )pbdoc",
         py::arg("cutoff") = 2.)
    .def("config_slipspring_low_cutoff",
         &dpd::DPDSimulator::configSlipspringLowCutoff,
         R"pbdoc(
          Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.
     )pbdoc",
         py::arg("cutoff") = 0.5)
    .def("config_box_deformation",
         &dpd::DPDSimulator::configBoxDeformation,
         R"pbdoc(
          Configure where to (incrementally) deform the box to during the next simulation run.
     )pbdoc",
         py::arg("target_box"))
#ifdef CEREALIZABLE
    .def_static("read_restart_file",
                &dpd::DPDSimulator::readRestartFile,
                R"pbdoc(
          Read a restart file in order to continue a simulation.
     )pbdoc",
                py::arg("file"))
    .def("config_restart_output",
         &dpd::DPDSimulator::configRestartOutput,
         R"pbdoc(
          Set when to output a restart where.

          Note:
               The filename determines the type of serialization:
               .json, .xml are supported; other file endings will lead to binary serialization (fastest!).

          Caution:
               This method may not be backwards- nor forward-compatible.
               Use the same version of pylimer-tools if you want to be sure that things work.

          Arguments:
               - file: the file path to the restart file to write
               - outputEvery: how often to write the restart file
     )pbdoc",
         py::arg("file"),
         py::arg("output_every") = 50000)
#endif
    .def("configAverageOutput",
         &dpd::DPDSimulator::configAverageOutput,
         R"pbdoc(
          Set which values to compute averages for.

          Arguments:
               - values: a list of OutputConfiguration structs
     )pbdoc")
    .def("configAutoCorrelatorOutput",
         &dpd::DPDSimulator::configAutoCorrelatorOutput,
         R"pbdoc(
          Set which values to compute multiple-tau autocorrelation for.
          If you use this, you should cite `doi:10.1063/1.3491098 <https://pubs.aip.org/aip/jcp/article-abstract/133/15/154103/190247/Efficient-on-the-fly-calculation-of-time?redirectedFrom=fulltext>`_

          Arguments:
               - values: a list of OutputConfiguration structs
               - ...
     )pbdoc",
         py::arg("values"),
         py::arg("num_corr_in") = 32,
         py::arg("p") = 16,
         py::arg("m") = 2)
    .def("config_step_output",
         &dpd::DPDSimulator::configStepOutput,
         R"pbdoc(
          Set which values to log.

          Arguments:
               - values: a list of OutputConfiguration structs
     )pbdoc")
    .def("config_shift_possibility_empty",
         &dpd::DPDSimulator::configShiftPossibilityEmpty,
         R"pbdoc()pbdoc")
    .def("config_shift_one_at_a_time",
         &dpd::DPDSimulator::configShiftOneAtATime,
         R"pbdoc()pbdoc")
    .def("config_num_steps_mc",
         &dpd::DPDSimulator::configNumStepsMC,
         R"pbdoc(
          Configure the number of steps to do in one MC sequence.
     )pbdoc")
    .def("config_num_steps_dpd",
         &dpd::DPDSimulator::configNumStepsDPD,
         R"pbdoc(
          Configure the number of steps to do in one DPD sequence.
     )pbdoc")
    .def("config_bond_formation",
         &dpd::DPDSimulator::configBondFormation,
         R"pbdoc(
          Configure how to do bond formation during the run.

          Arguments:
          - num_bonds_to_form (int): the nr of bonds to form in total. Use 0 to stop bond formation.
          - num_bonds_per_atom_type (dict): the nr of bonds each atom type may have at most (e.g., 2 for strand atoms, 4 for a tertiary crosslinkers)
          - bond_formation_dist (float): the maximum distance allowed to form bonds
          - attempt_bond_formation_every (int): attempt to form bonds every this many steps during the simulation run
          - atom_type_form_from (int): the atom type to start forming bonds from.
          - atom_type_form_to (int): the atom type to start forming bonds to.
         )pbdoc",
         py::arg("num_bonds_to_form"),
         py::arg("max_bonds_per_atom_type"),
         py::arg("bond_formation_dist") = 1.0,
         py::arg("attempt_bond_formation_every") = 50,
         py::arg("atom_type_form_from") = 2,
         py::arg("atom_type_form_to") = 1)
    .def("get_nr_of_bonds_to_form",
         &dpd::DPDSimulator::getNrOfBondsToForm,
         R"pbdoc(
          Get the number of bonds that are configured to have to be formed.
     )pbdoc")
    .def("config_allow_relocation_in_network",
         &dpd::DPDSimulator::configAllowRelocationInNetwork,
         R"pbdoc(
          Configure whether a relocation step may happen when a slip-spring has ended at a cross-link.

          Side-effect: if true, the relocations may also happen *to* a slip-spring next to a cross-link.

          Arguments:
          - allow_relocation_in_network (bool): Whether to allow relocation in the network or not.
         )pbdoc",
         py::arg("allow_relocation_in_network") = false)
    .def("start_measuring_msd_for_atoms",
         &dpd::DPDSimulator::startMeasuringMSDForAtoms,
         R"pbdoc(
          Set a new origin for measuing the mean square displacement for a specified set of atoms
         )pbdoc",
         py::arg("atom_ids"))
    .def("get_universe",
         &dpd::DPDSimulator::getUniverse,
         R"pbdoc(
     Get a universe instance from the current coordinates (and connectivity).

     Arguments:
          - with_slip_springs (bool): whether to include slip-springs in the returned universe.
    )pbdoc",
         py::arg("with_slipsprings") = true)
    .def("refresh_current_state",
         &dpd::DPDSimulator::refreshCurrentState,
         R"pbdoc(
          After re-configuring the force-field parameters,
          this method should be called to update the current stress tensor etc.
     )pbdoc")
    .def("get_timestep", &dpd::DPDSimulator::getTimestep)
    .def("get_current_timestep", &dpd::DPDSimulator::getCurrentTimestep)
    .def("get_temperature", &dpd::DPDSimulator::getTemperature)
    .def("get_bond_lengths", &dpd::DPDSimulator::getBondLengths)
    .def("get_coordinates", &dpd::DPDSimulator::getCoordinates)
    .def("get_spring_constant", &dpd::DPDSimulator::getSpringConstant)
    .def("get_shift_one_at_a_time", &dpd::DPDSimulator::getShiftOneAtATime)
    .def("get_nr_of_slip_springs", &dpd::DPDSimulator::getNumSlipSprings)
    .def("get_nr_of_atoms", &dpd::DPDSimulator::getNumAtoms)
    .def("get_nr_of_extra_atoms", &dpd::DPDSimulator::getNumExtraAtoms)
    .def("get_nr_of_bonds", &dpd::DPDSimulator::getNumBonds)
    .def("get_nr_of_extra_bonds", &dpd::DPDSimulator::getNumExtraBonds)
    .def("get_stress_tensor", &dpd::DPDSimulator::getStressTensor)
    .def("get_nr_of_steps_dpd", &dpd::DPDSimulator::getNumStepsDPD)
    .def("get_nr_of_steps_mc", &dpd::DPDSimulator::getNumStepsMC)
    .def("get_volume", &dpd::DPDSimulator::getVolume)
    .def("get_slip_spring_bond_type", &dpd::DPDSimulator::getSlipSpringBondType)
    .def("get_shift_possibility_empty",
         &dpd::DPDSimulator::getShiftPossibilityEmpty)
#ifdef CEREALIZABLE
    .def("write_restart_file",
         &dpd::DPDSimulator::writeRestartFile,
         R"pbdoc(
          Explicitly force the writing of a restart file, now!

          Arguments:
          - file (str): the file path and name of the restart file to be written.
               Can end in xml, json or anything else (-> binary).
     )pbdoc",
         py::arg("file"))
#endif
    .def("validate_neighbour_list", &dpd::DPDSimulator::validateNeighbourlist)
    .def("validate_state", &dpd::DPDSimulator::validateState);
}

#endif /* PYBIND_CALC_H */
