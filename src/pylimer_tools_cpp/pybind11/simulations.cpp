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
    PYBIND11_OVERRIDE_PURE(returntype,
                           /* Return type */
                           MEHPForceEvaluator,
                           /* Parent class */
                           evaluateForceSetGradient,
                           /* Name of function in C++
                              (must match Python name) */
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

  py::enum_<ComputedIntValues>(
    m, "ComputedIntValues", "Integer output quantities")
#define X(e, n)                                                                \
  .value(#e, ComputedIntValues::e, "Results in the output column " #n ".")
    COMPUTED_INT_VALUES
#undef X
    ;

  py::enum_<ComputedDoubleValues>(
    m, "ComputedDoubleValues", "Floating point output quantities")
#define X(e, n)                                                                \
  .value(#e, ComputedDoubleValues::e, "Results in the output column " #n ".")
    COMPUTED_DOUBLE_VALUES
#undef X
    ;

  py::class_<OutputConfiguration>(m,
                                  "OutputConfiguration",
                                  py::module_local(),
                                  R"pbdoc(
     A configuration object to configure the output values and frequency
     for simulation classes in this package.

     This class specifies which quantities to output and how often to write them
     during simulations.
    )pbdoc")
    .def(py::init<>(),
         R"pbdoc(
          Create a new OutputConfiguration instance.

          :return: A new OutputConfiguration object with default settings
         )pbdoc")
    .def_readwrite("int_values",
                   &OutputConfiguration::intValues,
                   R"pbdoc(
                    List of integer-valued quantities to output.
                    
                    Use :class:`~pylimer_tools_cpp.ComputedIntValues` enum to specify which integer quantities
                    should be computed and written to output.
                   )pbdoc")
    .def_readwrite("double_values",
                   &OutputConfiguration::doubleValues,
                   R"pbdoc(
                    List of double-valued quantities to output.
                    
                    Use :class:`~pylimer_tools_cpp.ComputedDoubleValues` enum to specify which floating-point quantities
                    should be computed and written to output.
                   )pbdoc")
    .def_readwrite("use_every",
                   &OutputConfiguration::useEvery,
                   R"pbdoc(
     For autocorrelation and averaging, how often to include values.

     Use a value of 1 to take average of or autocorrelate, respectively,
     all values encountered during the simulation or optimization procedure.
     )pbdoc")
    .def_readwrite("append",
                   &OutputConfiguration::append,
                   R"pbdoc(
     Whether to append to the file or truncate it
     )pbdoc")
    .def_readwrite("filename",
                   &OutputConfiguration::filename,
                   R"pbdoc(
      The path and name of the file to write to.
      An empty string ("") means standard output (console).
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

  py::enum_<mehp::ExitReason>(m,
                              "ExitReason",
                              R"pbdoc(
An enum representing the reason for exiting
the simulation or optimization procedure.)pbdoc")
#define X(e, n) .value(#e, mehp::ExitReason::e, "Exit reason: " #n ".")
    EXIT_REASONS
#undef X
    ;

  m.def("inverse_langevin",
        &mehp::langevin_inv,
        R"pbdoc(
     A somewhat accurate (for :math:`x \in (-1, 1)`) implementation of the inverse Langevin function.

     Source: https://scicomp.stackexchange.com/a/30251

     :param x: Input value in the range (-1, 1)
     :return: Inverse Langevin function value
  )pbdoc",
        py::arg("x"));

  ////////////////////////////////////////////////////////////////
  // MARK: Network structures
  py::class_<mehp::Network>(m,
                            "SimplifiedNetwork",
                            R"pbdoc(
     A more efficient structure of the network for use in MEHP,
     namely :obj:`~pylimer_tools_cpp.MEHPForceRelaxation`.
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
    .def_readonly("assume_box_large_enough",
                  &mehp::Network::assumeBoxLargeEnough)
    .def_readonly("assume_complete", &mehp::Network::assumeComplete);

  py::class_<mehp::ForceBalanceNetwork>(m,
                                        "SimplifiedBalanceNetwork",
                                        R"pbdoc(
     A more efficient structure of the network for use in MEHP force balance,
     namely :obj:`~pylimer_tools_cpp.MEHPForceBalance`, though also passable to
     namely :obj:`~pylimer_tools_cpp.MEHPForceBalance2`.
     Consists usually only of the cross- and slip-links (and their connectivity),
     i.e., no "normal strand beads" in between, in order to reduce the degrees of freedom
     and therewith improve performance of the solver.

     A note on the terminology: a spring is the connection between two links (crosslink, entanglement-link/slip-link).
     A strand is a chain of connected links between two crosslinks.
 )pbdoc")
    .def_readonly("box_lengths", &mehp::ForceBalanceNetwork::L)
    .def_readonly("volume", &mehp::ForceBalanceNetwork::vol)
    .def_readonly("nr_of_crosslinks", &mehp::ForceBalanceNetwork::nrOfNodes)
    .def_readonly("nr_of_links", &mehp::ForceBalanceNetwork::nrOfLinks)
    .def_readonly("nr_of_strands", &mehp::ForceBalanceNetwork::nrOfSprings)
    .def_readonly("nr_of_springs",
                  &mehp::ForceBalanceNetwork::nrOfPartialSprings)
    // .def_readonly("nrOfLoops", &mehp::Network::nrOfLoops)
    .def_readonly("coordinates", &mehp::ForceBalanceNetwork::coordinates)
    .def_readonly("old_atom_ids", &mehp::ForceBalanceNetwork::oldAtomIds)
    .def_readonly("strand_coordinate_index_a",
                  &mehp::ForceBalanceNetwork::springCoordinateIndexA)
    .def_readonly("strand_coordinate_index_b",
                  &mehp::ForceBalanceNetwork::springCoordinateIndexB)
    .def_readonly("spring_coordinate_index_a",
                  &mehp::ForceBalanceNetwork::springPartCoordinateIndexA)
    .def_readonly("spring_coordinate_index_b",
                  &mehp::ForceBalanceNetwork::springPartCoordinateIndexB)
    .def_readonly("strand_index_a", &mehp::ForceBalanceNetwork::springIndexA)
    .def_readonly("strand_index_b", &mehp::ForceBalanceNetwork::springIndexB)
    .def_readonly("spring_index_a",
                  &mehp::ForceBalanceNetwork::springPartIndexA)
    .def_readonly("spring_index_b",
                  &mehp::ForceBalanceNetwork::springPartIndexB)
    .def_readonly("link_is_sliplink",
                  &mehp::ForceBalanceNetwork::linkIsSliplink)
    .def_readonly("springs_of_strand",
                  &mehp::ForceBalanceNetwork::localToGlobalSpringIndex)
    .def_readonly("strands_of_link",
                  &mehp::ForceBalanceNetwork::springIndicesOfLinks)
    .def_readonly("links_of_strand",
                  &mehp::ForceBalanceNetwork::linkIndicesOfSprings)
    .def_readonly("nr_of_crosslink_swaps_endured",
                  &mehp::ForceBalanceNetwork::nrOfCrosslinkSwapsEndured)
    .def_readonly("strand_contour_length",
                  &mehp::ForceBalanceNetwork::springsContourLength)
    .def_readonly("strand_of_spring",
                  &mehp::ForceBalanceNetwork::partialToFullSpringIndex)
    .def_readonly("spring_box_offset",
                  &mehp::ForceBalanceNetwork::springPartBoxOffset)
    // .def_readonly("springIsActive", &mehp::Network::springIsActive)
    ;

  py::class_<mehp::ForceBalance2Network>(m,
                                         "SimplifiedBalance2Network",
                                         R"pbdoc(
A more efficient structure of the network for use in
:obj:`~pylimer_tools_cpp.MEHPForceBalance2`.
Consists usually only of the cross- and slip-links (and their connectivity),
i.e., no "normal strand beads" in between, in order to reduce the degrees of freedom
and therewith improve performance of the solver.

A note on the terminology: a spring is the connection between two links (crosslink, entanglement-link/slip-link).
A strand is a chain of connected links between two crosslinks.
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
    .def_readonly("link_is_entanglement",
                  &mehp::ForceBalance2Network::linkIsEntanglement)
    .def_readonly("spring_is_entanglement",
                  &mehp::ForceBalance2Network::springIsEntanglement)
    .def_readonly("spring_coordinate_index_a",
                  &mehp::ForceBalance2Network::springCoordinateIndexA)
    .def_readonly("spring_coordinate_index_b",
                  &mehp::ForceBalance2Network::springCoordinateIndexB)
    .def_readonly("spring_index_a", &mehp::ForceBalance2Network::springIndexA)
    .def_readonly("spring_index_b", &mehp::ForceBalance2Network::springIndexB)
    .def_readonly("spring_contour_length",
                  &mehp::ForceBalance2Network::springContourLength)
    .def_readonly("springs_of_strand",
                  &mehp::ForceBalance2Network::springIndicesOfStrand)
    .def_readonly("strands_of_link",
                  &mehp::ForceBalance2Network::strandIndicesOfLink)
    .def_readonly("links_of_strand",
                  &mehp::ForceBalance2Network::linkIndicesOfStrand)
    .def_readonly("strand_of_spring",
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

          :param springDistances: The three coordinate differences for one spring.
          :param i: The row index of the stress tensor
          :param j: The column index of the stress tensor
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
  //      :meth:`~pylimer_tools_cpp.MEHPForceEvaluator.evaluate_stress_contribution`.

  //      :param n: The dimensionality of the problem (the nr. of spring
  //      coordinates) :param springDistances: The sequential (x, y, z) spring
  //      distances :param displacements: The displacements from the original
  //      coordinates
  //           (accessible by
  //           :func:`~pylimer_tools_cpp.CustomMEHPForceEvaluator.getNetwork().coordinates`)
  //      :param gradientNeeded: whether the gradient should be computed and
  //      returned

  //      Returns:
  //           - force: The result of the force computation.
  //           - gradient: The result of the gradient computation.
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
     where :math:`r` is the spring [between crosslinkers] length.

     Recommended optimization algorithm: "LD_LBFGS"

     :param kappa: The spring constant :math:`\kappa`
    )pbdoc")
    .def(py::init<double>(), py::arg("kappa") = 1.0);

  py::class_<mehp::NonGaussianSpringForceEvaluator, mehp::MEHPForceEvaluator>(
    m,
    "NonGaussianSpringForceEvaluator",
    R"pbdoc(
     This is equal to a spring evaluator for Langevin chains.

     The force for a certain spring is given by:
     :math:`f = 0.5 \cdot \\frac{1}{l} \scriptL^{-1}(\frac{r}{N\cdot l})`,
     where :math:`r` is the spring [between crosslinkers] length
     and :math:`\scriptL^{-1}` the inverse langevin function.

     Please note that the inverse langevin is only approximated.

     Recommended optimization algorithm: "LD_MMA"

     :param kappa: The spring constant :math:`\kappa`
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
    A small simulation tool for quickly minimizing the force between the crosslinker beads.

    This is the first of three force relaxation methods available in this library.
    The relevant feature of this implementation is the configurable spring potential.
    Consequently, it offers a variety of configurable non-linear solvers using NLoptLib.

    Please cite :cite:t:`gusev_numerical_2019` if you use this method in your work.
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

          :param universe: The universe to simulate with
          :param crosslinker_type: The atom type of the crosslinkers. Needed to reduce the network.
          :param is2d: Whether to ignore the z direction.
          :param force_evaluator: The force evaluator to use
          :param kappa: The spring constant
          :param remove_2functional_crosslinkers: Whether to replace two-functional crosslinkers with a "normal" chain bead
          :param remove_dangling_chains: Whether to remove dangling chains before running the simulation.
               **Caution**: Removing the dangling chains will result in incorrect results fo the computation of
               :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_soluble_weight_fraction` and
               :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_dangling_weight_fraction`
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
          
          :param force_evaluator: The new force evaluator to use
     )pbdoc",
         py::arg("force_evaluator"))
    .def("config_rerun_epsilon",
         &mehp::MEHPForceRelaxation::configRerunEps,
         R"pbdoc(
          Configure the offset from the lower and upper bounds for the simulation to suggest another run.
          
          :param epsilon: The epsilon value to use for the rerun check
               (See: :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.requires_another_run`)
         )pbdoc",
         py::arg("epsilon") = 1e-3)
    .def("config_step_output",
         &mehp::MEHPForceRelaxation::configStepOutput,
         R"pbdoc(
          Set which values to log during the simulation.

          :param output_configuration: An OutputConfiguration struct or list of OutputConfiguration structs
                                      specifying what values to log and how often
     )pbdoc",
         py::arg("output_configuration"))
    .def("assume_box_large_enough",
         &mehp::MEHPForceRelaxation::configAssumeBoxLargeEnough,
         R"pbdoc(
          Configure whether to run PBC on the bonds or not.

          :param box_large_enough: If True, assume the box is large enough and don't apply PBC on bonds.
                                  If your bonds could get larger than half the box length, 
                                  this must be kept False (default).
         )pbdoc",
         py::arg("box_large_enough") = false)
    .def("get_force",
         &mehp::MEHPForceRelaxation::getForce,
         R"pbdoc(
          Returns the force at the current state of the simulation.
          
          :return: The current force value
     )pbdoc")
    .def("get_residuals",
         &mehp::MEHPForceRelaxation::getResiduals,
         R"pbdoc(
          Returns the residuals at the current state of the simulation.

          :return: The current residual vector
     )pbdoc")
    .def("get_residual_norm",
         &mehp::MEHPForceRelaxation::getResidualNorm,
         R"pbdoc(
          Returns the residual norm at the current state of the simulation.

          :return: The current residual norm value
     )pbdoc")
    .def("get_pressure",
         &mehp::MEHPForceRelaxation::getPressure,
         R"pbdoc(
          Returns the pressure at the current state of the simulation.

          :return: The current pressure value
     )pbdoc")
    .def("get_stress_tensor",
         &mehp::MEHPForceRelaxation::getStressTensor,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.

          :return: The current stress tensor matrix
     )pbdoc")
    .def("get_gamma_factors",
         &mehp::MEHPForceRelaxation::getGammaFactors,
         R"pbdoc(
          Computes the gamma factor for each spring as part of the ANT/MEHP formulism.

          :math:`\gamma_{\eta} = \\frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`, with (here)
          :math:`R_{0,\eta}^2 = N_\eta \cdot ` the parameter `b0_squared`.
          You can obtain this parameter e.g. by doing melt simulations at different lengths,
          it's the slope you obtain.

          :param b0_squared: Part of the denominator in the equation of :math:`\Gamma`.
               If :math:`-1.0` (default), the network is used for determination (which is not accurate), the system is assumed to be phantom.
               For real systems, the value could be determined by :func:`~pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance()`
               on the melt system, with subsequent division by the nr of bonds in the chain.

          See also :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor` for the mean of these.
         )pbdoc",
         py::arg("b0_squared") = -1.0)
    .def("get_gamma_factor",
         &mehp::MEHPForceRelaxation::getGammaFactor,
         R"pbdoc(
          Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

          :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \\frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
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
          :param nr_of_chains: The value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains.
     )pbdoc",
         py::arg("b0_squared") = -1.0,
         py::arg("nr_of_chains") = -1)
    .def("get_nr_of_nodes",
         &mehp::MEHPForceRelaxation::getNrOfNodes,
         R"pbdoc(
          Get the number of nodes considered in this simulation.

          :return: The number of nodes in the simulation
     )pbdoc")
    .def("get_nr_of_springs",
         &mehp::MEHPForceRelaxation::getNrOfSprings,
         R"pbdoc(
          Get the number of springs considered in this simulation.

          :param tolerance: Springs under this length are considered inactive
          :return: The number of springs in the simulation
     )pbdoc")
    .def("get_ids_of_active_nodes",
         &mehp::MEHPForceRelaxation::getIdsOfActiveNodes,
         R"pbdoc(
          Get the atom ids of the nodes that are considered active.

          :param tolerance: Springs under this length are considered inactive. A node is active if it has > 2 active springs.
          :param minimum_nr_of_active_connections: Minimum number of active connections required for a node to be considered active
          :param maximum_nr_of_active_connections: Maximum number of active connections allowed for a node to be considered active
          :return: List of atom IDs for active nodes
     )pbdoc",
         py::arg("tolerance") = 1e-3,
         py::arg("minimum_nr_of_active_connections") = 2,
         py::arg("maximum_nr_of_active_connections") = -1)
    .def("get_nr_of_active_nodes",
         &mehp::MEHPForceRelaxation::getNrOfActiveNodes,
         R"pbdoc(
          Get the number of active nodes remaining after running the simulation.

          :param tolerance: Springs under this length are considered inactive
          :param minimumNrOfActiveConnections: A node is active if it has equal or more than this number of active springs
          :param maximumNrOfActiveConnections: A node is active if it has equal or less than this number of active springs.
               Use a value < 0 to indicate that there is no maximum number of active connections
          :return: The number of active nodes
     )pbdoc",
         py::arg("tolerance") = 1e-3,
         py::arg("minimumNrOfActiveConnections") = 2,
         py::arg("maximumNrOfActiveConnections") = -1)
    .def("get_nr_of_active_springs",
         &mehp::MEHPForceRelaxation::getNrOfActiveSprings,
         R"pbdoc(
          Get the number of active springs remaining after running the simulation.

          :param tolerance: Springs under this length are considered inactive
          :return: The number of active springs
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("count_active_clustered_atoms",
         py::overload_cast<const double>(
           &mehp::MEHPForceRelaxation::countActiveClusteredAtoms),
         R"pbdoc(
          Counts the active clustered atoms in the system.

          :param tolerance: Springs under this length are considered inactive
          :return: The number of active clustered atoms
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_soluble_weight_fraction",
         &mehp::MEHPForceRelaxation::getSolubleWeightFraction,
         R"pbdoc(
          Compute the weight fraction of springs connected to active
          springs (any depth).

          Caution: ignores atom masses.

          :param tolerance: Springs under this length are considered inactive
          :return: The soluble weight fraction
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_dangling_weight_fraction",
         &mehp::MEHPForceRelaxation::getDanglingWeightFraction,
         R"pbdoc(
          Compute the weight fraction of non-active springs.

          Caution: ignores atom masses.

          :param tolerance: Springs under this length are considered inactive
          :return: The dangling weight fraction
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_active_chains",
         &mehp::MEHPForceRelaxation::getActiveChains,
         R"pbdoc(
          Get the crosslinker chains that are active.

          :param tolerance: Springs under this length are considered inactive
          :return: List of active crosslinker chains
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_effective_functionality_of_atoms",
         &mehp::MEHPForceRelaxation::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
          Returns the number of active springs connected to each atom, atomId used as index.

          :param tolerance: Springs under this length are considered inactive
          :return: Vector with effective functionality for each atom (indexed by atom ID)
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_spring_lengths",
         &mehp::MEHPForceRelaxation::getSpringLengths,
         R"pbdoc(
          Get the current lengths for all the springs.

          :return: A vector of size nrOfSprings, with each the norm of the distances
     )pbdoc")
    .def("get_spring_distances",
         &mehp::MEHPForceRelaxation::getSpringDistances,
         R"pbdoc(
          Get the current coordinate differences for all the springs.

          :return: A vector of size 3*nrOfSprings, with each x, y, z values of the springs
     )pbdoc")
    .def("get_average_spring_length",
         &mehp::MEHPForceRelaxation::getAverageSpringLength,
         R"pbdoc(
          Get the average length of the springs. Note that in contrast to :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor`,
          this value is normalized by the number of springs rather than the number of chains.

          :return: The average spring length
     )pbdoc")
    .def("get_default_r0_square",
         &mehp::MEHPForceRelaxation::getDefaultR0Square,
         R"pbdoc(
          Returns the value effectively used in :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor` for :math:`\langle R_{0,\eta}^2\rangle`.

          :return: The default R0 squared value used in gamma factor calculations
     )pbdoc")
    .def("get_nr_of_iterations",
         &mehp::MEHPForceRelaxation::getNrOfIterations,
         R"pbdoc(
          Returns the number of iterations used for force relaxation.

          :return: The number of iterations performed during force relaxation
     )pbdoc")
    .def("get_exit_reason",
         &mehp::MEHPForceRelaxation::getExitReason,
         R"pbdoc(
          Returns the reason for termination of the simulation.

          :return: The exit reason enum value indicating why the simulation ended
     )pbdoc")
    .def("requires_another_run",
         &mehp::MEHPForceRelaxation::suggestsRerun,
         R"pbdoc(
          For performance reasons, the objective is only minimised within the distances of one box.
          This means, that there is a possibility, e.g. for a single strand longer than two boxes,
          that it would not be globally minimised.

          If the final displacement of one of the atoms is close
          (1e-3, configurable via :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.config_rerun_epsilon`)
          to the imposed min/max, after minimizing,
          this method would return true.

          :return: True if another run is suggested, False otherwise
     )pbdoc")
    .def("get_crosslinker_universe",
         &mehp::MEHPForceRelaxation::getCrosslinkerVerse,
         R"pbdoc(
          Returns the universe [of crosslinkers] with the positions of the current state of the simulation.

          :return: A Universe object containing the crosslinkers with updated positions
     )pbdoc")
    .def_property_readonly("network",
                           &mehp::MEHPForceRelaxation::getNetwork,
                           R"pbdoc(
          Get the network structure.

          :return: The network structure used in the simulation
     )pbdoc")
    .def("get_spring_contour_length",
         &mehp::MEHPForceRelaxation::getSpringContourLength,
         R"pbdoc(
          Get the spring contour lengths.

          :return: A vector containing the contour lengths of all springs
     )pbdoc")
    .def("get_nr_of_active_springs_connected",
         &mehp::MEHPForceRelaxation::getNrOfActiveSpringsConnected,
         R"pbdoc(
          Returns the number of active springs connected to each node.

          :param tolerance: Springs under this length are considered inactive
          :return: Vector with the number of active springs for each node
     )pbdoc",
         py::arg("tolerance") = 0.05)
    .def("get_current_spring_distances",
         &mehp::MEHPForceRelaxation::getCurrentSpringDistances,
         R"pbdoc(
          Get the current spring distances.

          :return: A vector containing the current spring distance vectors
     )pbdoc")
    .def("get_average_contour_length",
         &mehp::MEHPForceRelaxation::getAverageContourLength,
         R"pbdoc(
          Get the average contour length of all springs.

          :return: The average contour length
     )pbdoc")
    .def("config_step_output",
         &mehp::MEHPForceRelaxation::configStepOutput,
         R"pbdoc(
          Set which values to log during the simulation.

          :param output_configuration: An OutputConfiguration struct or list of OutputConfiguration structs
                                      specifying what values to log and how often
     )pbdoc",
         py::arg("output_configuration"))
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
  py::enum_<mehp::StructureSimplificationMode>(
    m,
    "StructureSimplificationMode",
    "How the structure shall be simplified during the optimization in order to "
    "remove non-trapped entanglement links.")
#define X(e, name)                                                             \
  .value(STRINGINFY(e), mehp::StructureSimplificationMode::e, name)
    STRUCTURE_SIMPLIFICATION_MODES;
#undef X

  py::enum_<mehp::LinkSwappingMode>(
    m,
    "LinkSwappingMode",
    "How slip-links may act when they reach each-other or even a crosslink.")
#define X(e, name) .value(STRINGINFY(e), mehp::LinkSwappingMode::e, name)
    LINK_SWAPPING_MODES;
#undef X

  py::enum_<mehp::SLESolver>(
    m, "SLESolver", "Solver for sparse linear equation systems")
#define X(e, name) .value(STRINGINFY(e), mehp::SLESolver::e, name)
    SLE_SOLVERS;
#undef X

  ////////////////////////////////////////////////////////////////
  // MARK: Force Balance
  py::class_<mehp::MEHPForceBalance>(m,
                                     "MEHPForceBalance",
                                     R"pbdoc(
    A small simulation tool for quickly minimizing the force between the crosslinker beads.

    This is the second implementation in the group of MEHP provided by this package.
    The distinct feature here is the slip-links: a form of entanglement,
    represented as an entanglement link, just like a four-functional crosslink,
    but with the ability to slip along the two associated strands, therewith
    adjusting the fraction of the contour length on both sides of the link.

    Please cite :cite:t:`bernhard_phantom_2025` if you use this method.
    )pbdoc")
    .def(py::init<pe::Universe, int, bool, bool, bool>(),
         R"pbdoc(
          Instantiate the simulator for a certain universe.

          :param universe: The universe to simulate with
          :param crosslinker_type: The atom type of the crosslinkers. Needed to reduce the network.
          :param is2D: Whether to ignore the z direction.
          :param kappa: The spring constant
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

          :return: True if validation passes, raises exception otherwise
     )pbdoc")
    .def(
      "run_force_relaxation",
      [](mehp::MEHPForceBalance& sim,
         const long int maxNrOfSteps,
         // default: 10000
         const double xtol,
         const double initialResidualToUse,
         const mehp::StructureSimplificationMode simplificationMode,
         const double inactiveRemovalCutoff,
         const bool doInnerIterations,
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
          :param one_over_strand_partition_upper_limit: Super-secret parameter. Use 1, gradually increase (and then -1) if you want to publish.
          :param nr_of_crosslink_swaps_allowed_per_sliplink: Use to steer whether slip-links can cross crosslinks when swapping is enabled.
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
      py::arg("one_over_strand_partition_upper_limit") = 1.0,
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
          Set which values to log during the simulation.

          :param output_configuration: A list of OutputConfiguration structs
                                      specifying what values to log and how often
     )pbdoc",
         py::arg("output_configuration"))
    .def("config_assume_box_large_enough",
         &mehp::MEHPForceBalance::configAssumeBoxLargeEnough,
         R"pbdoc(
          Configure whether to run PBC on the bonds or not.

          :param box_large_enough: If True, assume the box is large enough and don't apply PBC on bonds.
                                  If your bonds could get larger than half the box length, 
                                  this must be kept False (default).
         )pbdoc",
         py::arg("box_large_enough") = false)
    .def("config_mean_bond_length",
         &mehp::MEHPForceBalance::configMeanBondLength,
         R"pbdoc(
          Configure the :math:`b` used e.g. for the topological Gamma-factor.

          :param b: The mean bond length to use.
     )pbdoc",
         py::arg("b") = 1.0)
    .def("config_spring_constant",
         &mehp::MEHPForceBalance::configSpringConstant,
         R"pbdoc(
          Configure the spring constant used in the simulation.

          :param kappa: The spring constant.
     )pbdoc",
         py::arg("kappa") = 1.0)
    .def("config_entanglement_type",
         &mehp::MEHPForceBalance::configEntanglementType,
         R"pbdoc(
         To have certain crosslinks behave as entanglements in the removal process,
         you can specify here a type, that you have used in the universe to specify:
         - the type of entanglement atoms (expected with functionality f = 3),
         - and the entanglement-bonds between the entanglement atoms.

         I.e., say you want to model some entanglements as non-slipping,
         bonds between two strand beads resulting in f = 3 beads, for example,
         you can call this method to have the "StructureSimplificationMode" also remove these atoms,
         if they have a functionality of 2 or less while still being connected to its partner bead.

         :param type: The atom type to treat as entanglement.
     )pbdoc",
         py::arg("type") = -1)
    .def("config_spring_breaking_distance",
         &mehp::MEHPForceBalance::configSpringBreakingDistance,
         R"pbdoc(
          Configure the "force" (distance over contour length) at which the bonds break.
          Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
          Springs breaking will happen before the simplification procedure is run.
          Negative values will disable spring breaking.
          Default: -1.

          :param distance_over_contour_length: The threshold for breaking springs.
         )pbdoc",
         py::arg("distance_over_contour_length") = -1)
    .def("config_simplification_frequency",
         &mehp::MEHPForceBalance::configSimplificationFrequency,
         R"pbdoc(
          Configure every how many steps to simplify the structure.

          :param frequency: The number of steps between simplification. Default: 10.
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
      py::arg("one_over_strand_partition_upper_limit") = 1.0)
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
      py::arg("one_over_strand_partition_upper_limit") = 1.0)
    .def("inspect_displacement_to_mean_position_update",
         &mehp::MEHPForceBalance::inspectDisplacementToMeanPositionUpdate,
         R"pbdoc(
          Helper method to debug and/or understand what happens to certain links when being displaced.
         )pbdoc",
         py::arg("link_idx"),
         py::arg("one_over_strand_partition_upper_limit") = 1.0)
    .def("inspect_strand_partition_update",
         &mehp::MEHPForceBalance::inspectSpringPartitionUpdate,
         R"pbdoc(
          Helper method to debug and/or understand what happens to certain links
          when the strand partition is being updated.
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
         py::arg("strand_partitions"),
         py::arg("max_nr_of_steps") = 100,
         py::arg("alpha_tol") = 1e-9,
         py::arg("min_nr_of_steps") = 1,
         py::arg("one_over_strand_partition_upper_limit") = 1.0)
    .def("get_strand_partition_indices_of_sliplink",
         &mehp::MEHPForceBalance::getSpringpartitionIndicesOfSliplink,
         R"pbdoc(
          Get the indices of strand partitions associated with a slip-link.

          :param network: The force balance network
          :param link_idx: Index of the slip-link
          :return: Vector of strand partition indices
         )pbdoc",
         py::arg("network"),
         py::arg("link_idx"))
    .def_static("get_neighbour_link_indices",
                &mehp::MEHPForceBalance::getNeighbourLinkIndices,
                R"pbdoc(
                 Get the indices of neighboring links for a given link.

                 :param network: The force balance network
                 :param link_idx: Index of the link to find neighbors for
                 :return: Vector of connected link indices
                )pbdoc",
                py::arg("network"),
                py::arg("link_idx"))
    .def(
      "evaluate_spring_distance",
      [](const mehp::MEHPForceBalance& sim,
         const mehp::ForceBalanceNetwork& net,
         const Eigen::VectorXd& u,
         const size_t springIdx) {
        return sim.evaluatePartialSpringDistance(net, u, springIdx);
      },
      R"pbdoc(
       Evaluate the distance vector for a specific spring.

       :param network: The force balance network
       :param displacements: Current displacement vector
       :param spring_idx: Index of the spring to evaluate
       :return: 3D distance vector for the spring
      )pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"))
    .def(
      "evaluate_spring_distance_from",
      [](const mehp::MEHPForceBalance& sim,
         const mehp::ForceBalanceNetwork& net,
         const Eigen::VectorXd& u,
         const size_t springIdx,
         const size_t linkIdx) {
        return sim.evaluatePartialSpringDistanceFrom(
          net, u, springIdx, linkIdx);
      },
      R"pbdoc(
       Evaluate the spring distance from a specific link.

       :param network: The force balance network
       :param displacements: Current displacement vector
       :param spring_idx: Index of the spring to evaluate
       :param link_idx: Index of the starting link
       :return: 3D distance vector from the specified link
      )pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"),
      py::arg("link_idx"))
    .def(
      "evaluate_spring_distance_to",
      [](const mehp::MEHPForceBalance& sim,
         const mehp::ForceBalanceNetwork& net,
         const Eigen::VectorXd& u,
         const size_t springIdx,
         const size_t linkIdx) {
        return sim.evaluatePartialSpringDistanceTo(net, u, springIdx, linkIdx);
      },
      R"pbdoc(
       Evaluate the spring distance to a specific link.

       :param network: The force balance network
       :param displacements: Current displacement vector
       :param spring_idx: Index of the spring to evaluate
       :param link_idx: Index of the target link
       :return: 3D distance vector to the specified link
      )pbdoc",
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
          Compute the weight fraction of strands connected to active
          strands (any depth).

          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_dangling_weight_fraction",
         &mehp::MEHPForceBalance::getDanglingWeightFraction,
         R"pbdoc(
          Compute the weight fraction of non-active strands

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
      py::arg("one_over_strand_partition_upper_limit") = 1.)
    .def("get_stress_tensor_link_based",
         &mehp::MEHPForceBalance::getStressTensorLinkBased,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc",
         py::arg("one_over_strand_partition_upper_limit") = 1.,
         py::arg("xlinks_only") = false)
    .def("get_gamma_factor",
         &mehp::MEHPForceBalance::getGammaFactor,
         R"pbdoc(
          Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

          :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \\frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
          which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
          where :math:`\eta` is the index of a particular strand,
          :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\eta}*b^2$`
          :math:`N_{\eta}` is the number of atoms in this strand :math:`\eta`,
          :math:`b` its mean square bond length,
          :math:`\nu` the number density of network strands,
          :math:`T` the temperature and
          :math:`k_B` Boltzmann's constant.

          :param b02: The melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\infinity b^2`.
          :param nr_of_chains: The value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of springs.
     )pbdoc",
         py::arg("b02") = -1.0,
         py::arg("nr_of_chains") = -1,
         py::arg("one_over_strand_partition_upper_limit") = 1.)
    .def("get_gamma_factors",
         &mehp::MEHPForceBalance::getGammaFactors,
         R"pbdoc(
          Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)
     )pbdoc",
         py::arg("b02"),
         py::arg("one_over_strand_partition_upper_limit") = 1.)
    .def("get_gamma_factors_in_dir",
         &mehp::MEHPForceBalance::getGammaFactorsInDir,
         R"pbdoc(
               Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)

               :param b02: The melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\infinity b^2`.
               :param direction: The direction in which to compute the gamma factors (0: x, 1: y, 2: z)
          )pbdoc",
         py::arg("b02"),
         py::arg("direction"),
         py::arg("one_over_strand_partition_upper_limit") = 1.)
    .def("get_nr_of_nodes",
         &mehp::MEHPForceBalance::getNrOfNodes,
         R"pbdoc(
           Get the number of nodes (crosslinkers) considered in this simulation.
     )pbdoc")
    .def("get_nr_of_strands",
         &mehp::MEHPForceBalance::getNrOfSprings,
         R"pbdoc(
          Get the number of strands considered in this simulation.
     )pbdoc")
    .def("get_strand_partitions",
         &mehp::MEHPForceBalance::getSpringPartitions,
         R"pbdoc(
          Get the current strand partitions (the fraction of the contour length associated with each spring).
     )pbdoc")
    .def("get_weighted_spring_lengths",
         &mehp::MEHPForceBalance::getWeightedPartialSpringLengths,
         R"pbdoc(
          Get the current spring lengths (norm of vector) divided by the strand partition times the contour length.
          )pbdoc",
         py::arg("one_over_strand_partition_upper_limit") = 1.)
    .def("set_strand_partitions",
         &mehp::MEHPForceBalance::setSpringPartitions,
         R"pbdoc(
          Set the current strand partitions.
     )pbdoc")
    .def("get_coordinates",
         &mehp::MEHPForceBalance::getCoordinates,
         R"pbdoc(
             Get the current coordinates of the crosslinkers and entanglement links.
      )pbdoc")
    .def("get_initial_coordinates",
         &mehp::MEHPForceBalance::getInitialCoordinates,
         R"pbdoc(
             Get the initial coordinates of the (remaining) crosslinkers and entanglement links.
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
    .def("set_strand_contour_lengths",
         &mehp::MEHPForceBalance::setSpringContourLengths,
         R"pbdoc(
           Set/overwrite the contour lengths.
      )pbdoc")
    .def("get_displacement_residual_norm",
         &mehp::MEHPForceBalance::getDisplacementResidualNorm,
         R"pbdoc(
           Get the current link displacement residual norm.
      )pbdoc",
         py::arg("one_over_strand_partition_upper_limit") = 1.)
    .def("get_ids_of_active_nodes",
         &mehp::MEHPForceBalance::getIdsOfActiveNodes,
         R"pbdoc(
          Get the atom ids of the nodes that are considered active.
          Only crosslink ids are returned (not e.g. entanglement links).

          :param tolerance: strands under this length are considered inactive. A node is active if it has > 1 active strands.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_nodes",
         &mehp::MEHPForceBalance::getNrOfActiveNodes,
         R"pbdoc(
          Get the number of active nodes (incl. entanglement nodes [atoms with type = entanglementType, present in the universe when creating this simulator],
          excl. entanglement links [the slip-links created internally when e.g. constructing the simulator with random slip-links]).

          :param tolerance: strands under this length are considered inactive. A node is active if it has > 1 active strands.
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_strands",
         &mehp::MEHPForceBalance::getNrOfActiveSprings,
         R"pbdoc(
           Get the number of active strands remaining after running the simulation.

          :param tolerance: strands under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_strands_in_dir",
         &mehp::MEHPForceBalance::getNrOfActiveSpringsInDir,
         R"pbdoc(
                Get the number of active strands remaining after running the simulation.

               :param direction: The direction in which to compute the active strands (0: x, 1: y, 2: z)
               :param tolerance: strands under this length are considered inactive
          )pbdoc",
         py::arg("direction"),
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_springs",
         &mehp::MEHPForceBalance::getNrOfActivePartialSprings,
         R"pbdoc(
           Get the number of active springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_current_spring_vectors",
         &mehp::MEHPForceBalance::getCurrentPartialSpringDistances,
         R"pbdoc(
          Get the spring vectors.
     )pbdoc")
    .def("get_current_spring_lengths",
         &mehp::MEHPForceBalance::getCurrentPartialSpringLengths,
         R"pbdoc(
          Get the spring distances.
     )pbdoc")
    .def("get_current_strand_vectors",
         &mehp::MEHPForceBalance::getCurrentSpringDistances)
    .def("get_overall_strand_lengths",
         &mehp::MEHPForceBalance::getOverallSpringLengths,
         R"pbdoc(
          Get the sum of the lengths of the springs of each strand.
     )pbdoc")
    .def("get_effective_functionality_of_atoms",
         &mehp::MEHPForceBalance::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
          Returns the number of active strands connected to each atom, atomId used as index

          :param tolerance: strands under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_average_strand_length",
         &mehp::MEHPForceBalance::getAverageSpringLength,
         R"pbdoc(
           Get the average length of the strands. Note that in contrast to :meth:`~pylimer_tools_cpp.MEHPForceBalance.get_gamma_factor`,
           this value is normalized by the number of strands rather than the number of chains.
     )pbdoc")
    .def("get_default_mean_bond_length",
         &mehp::MEHPForceBalance::getDefaultMeanBondLength,
         R"pbdoc(
           Returns the value effectively used in :meth:`~pylimer_tools_cpp.MEHPForceBalance.get_gamma_factor` for
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
           Returns the universe [of crosslinkers] with the positions of the current state of the simulation.
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
     A small simulation tool for quickly minimizing the force between the crosslinker beads.

     This is the third implementation of the MEHP. 
     It's the fastest implementation by using a simple spring model only, disabling the non-linear
     periodic boundary conditions, and instead builds a sparse linear system of equations that's readily solved.
     However, it allows entanglements to be represented as additional links or/and springs,
     although without slipping along the chain.

     Please cite :cite:t:`bernhard_phantom_2025` if you use this method.
      )pbdoc")
    .def(py::init<pe::Universe, int, bool>(),
         R"pbdoc(
           Instantiate the simulator for a certain universe.

           :param universe: The universe giving the basic connectivity to compute with
           :param crosslinker_type: The atom type of the crosslinkers. Needed to reduce the network.
           :param is_2d: Whether to ignore the z direction.
           :param kappa: The spring constant
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

           :param universe: The universe giving the basic connectivity to compute with
           :param entanglements: The entanglements to use in the computation
           :param crosslinker_type: The atom type of the crosslinkers. Needed to reduce the network.
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

     :param universe: The universe containing the basic atoms and connectivity
     :param nr_of_entanglements_to_sample: The number of entanglements to sample
     :param upper_cutoff: maximum distance from one sampled bead to its partner
     :param lower_cutoff: minimum distance from one sampled bead to its partner
     :param minimum_nr_of_entanglements: The minimum number of entanglements that should be sampled
     :param same_strand_cutoff: distance from one sampled bead to its pair within the same strand
     :param seed: The seed for the random number generator
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

           :return: True if validation passes, raises exception otherwise
      )pbdoc")
    .def(
      "run_force_relaxation",
      [](mehp::MEHPForceBalance2& sim,
         const mehp::StructureSimplificationMode simplificationMode,
         const double inactiveRemovalCutoff,
         const mehp::SLESolver solverChoice,
         const double tolerance,
         const int maxIterations) {
        return sim.runForceRelaxation(
          simplificationMode,
          inactiveRemovalCutoff,
          solverChoice,
          tolerance,
          maxIterations,
          []() { return PyErr_CheckSignals() != 0; },
          []() { throw py::error_already_set(); });
      },
      R"pbdoc(
           Run the simulation.

           :param simplification_mode: How to simplify the structure during the minimization.
           :param inactive_removal_cutoff: The tolerance in distance units of the partial spring length to count as active, relevant if simplification mode is specified to be something other than NO_SIMPLIFICATION.
           :param sle_solver: The solver to use for the system of linear equations (SLE).
           :param tolerance: The stopping condition/tolerance for the SLE solver if it's an iterative one.
           :param max_iterations: The maximum number of iterations to perform if the SLE solver is an iterative one.
           )pbdoc",
      py::arg("simplification_mode") =
        mehp::StructureSimplificationMode::NO_SIMPLIFICATION,
      py::arg("inactive_removal_cutoff") = 1e-6,
      py::arg("sle_solver") = mehp::SLESolver::DEFAULT,
      py::arg("tolerance") = 1e-15,
      py::arg("max_iterations") = 10000)
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

           :param values: a list of OutputConfiguration structs
      )pbdoc",
         py::arg("values"))
    .def("config_mean_bond_length",
         &mehp::MEHPForceBalance2::configMeanBondLength,
         R"pbdoc(
      Configure the :math:`b` used e.g. for the topological Gamma-factor.

      :param b: The mean bond length to use.
      )pbdoc",
         py::arg("b") = 1.0)
    .def("config_spring_constant",
         &mehp::MEHPForceBalance2::configSpringConstant,
         R"pbdoc(
           Configure the spring constant used in the simulation.

           :param kappa: The spring constant.
      )pbdoc",
         py::arg("kappa") = 1.0)
    .def("config_spring_breaking_distance",
         &mehp::MEHPForceBalance2::configSpringBreakingDistance,
         R"pbdoc(
           Configure the "force" (distance over contour length) at which the bonds break.
           Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
           Springs breaking will happen before the simplification procedure is run.
           Negative values will disable spring breaking.
           Default: -1.

           :param distance_over_contour_length: The threshold for breaking springs.
          )pbdoc",
         py::arg("distance_over_contour_length") = -1)
    .def("config_assume_network_is_complete",
         &mehp::MEHPForceBalance2::configAssumeNetworkIsComplete,
         R"pbdoc(
           Configure whether the network is assumed to be complete.

           This assumption means, that the `universe` instance is not queried for clusters when
           computing the fraction of e.g. the soluble or dangling strands.
           Do not set this to true if inactive
           (dangling or free) strands have been removed from the network.

           This is useful to reduce memory between MC generator and force balance,
           since the universe representation with all the in-between beads does not need to be stored.
           However, currently, the removal procedures are not loss-free.
           Therefore, use this with care.
         )pbdoc",
         py::arg("assume_network_is_complete") = false)
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
    .def_static("get_neighbour_link_indices",
                &mehp::MEHPForceBalance2::getNeighbourLinkIndices,
                R"pbdoc(
                 Get the indices of neighboring links for a given link.

                 :param network: The force balance 2 network
                 :param link_idx: Index of the link to find neighbors for
                 :return: Vector of connected link indices
                )pbdoc",
                py::arg("network"),
                py::arg("link_idx"))
    .def(
      "evaluate_spring_vector",
      [](const mehp::MEHPForceBalance2& sim,
         const mehp::ForceBalance2Network& net,
         const Eigen::VectorXd& u,
         const size_t springIdx) {
        return sim.evaluateSpringVector(net, u, springIdx);
      },
      R"pbdoc(
       Evaluate the distance vector for a specific spring.

       :param network: The force balance 2 network
       :param displacements: Current displacement vector
       :param spring_idx: Index of the spring to evaluate
       :return: 3D distance vector for the spring
      )pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"))
    .def(
      "evaluate_spring_vector_from",
      [](const mehp::MEHPForceBalance2& sim,
         const mehp::ForceBalance2Network& net,
         const Eigen::VectorXd& u,
         const size_t springIdx,
         const size_t linkIdx) {
        return sim.evaluateSpringVectorFrom(net, u, springIdx, linkIdx);
      },
      R"pbdoc(
       Evaluate the spring vector in the direction from a specific link.

       :param network: The force balance 2 network
       :param displacements: Current displacement vector
       :param spring_idx: Index of the spring to evaluate
       :param link_idx: Index of the starting link
       :return: 3D distance vector from the specified link
      )pbdoc",
      py::arg("network"),
      py::arg("displacements"),
      py::arg("spring_idx"),
      py::arg("link_idx"))
    .def(
      "evaluate_spring_vector_to",
      [](const mehp::MEHPForceBalance2& sim,
         const mehp::ForceBalance2Network& net,
         const Eigen::VectorXd& u,
         const size_t springIdx,
         const size_t linkIdx) {
        return sim.evaluateSpringVectorTo(net, u, springIdx, linkIdx);
      },
      R"pbdoc(
       Evaluate the spring vector in the direction to a specific link.

       :param network: The force balance 2 network
       :param displacements: Current displacement vector
       :param spring_idx: Index of the spring to evaluate
       :param link_idx: Index of the target link
       :return: 3D distance vector to the specified link
      )pbdoc",
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

           :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \\frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
           which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
           where :math:`\eta` is the index of a particular strand,
           :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\eta}*b^2$`
           :math:`N_{\eta}` is the number of atoms in this strand :math:`\eta`,
           :math:`b` its mean square bond length,
           :math:`T` the temperature and
           :math:`k_B` Boltzmann's constant.

           :param b02: The melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\infinity b^2`.
           :param nr_of_chains: The value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of springs.
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

                :param b02: The melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\infinity b^2`.
                :param direction: The direction in which to compute the gamma factors (0: x, 1: y, 2: z)
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
    .def("get_weighted_spring_lengths",
         &mehp::MEHPForceBalance2::getWeightedSpringLengths,
         R"pbdoc(
           Get the current spring lengths (norm of vector) divided by the spring partition times the contour length.
           )pbdoc")
    .def("get_initial_coordinates",
         &mehp::MEHPForceBalance2::getInitialCoordinates,
         R"pbdoc(
             Get the initial coordinates of the crosslinkers and entanglement links.
      )pbdoc")
    .def("get_coordinates",
         &mehp::MEHPForceBalance2::getCoordinates,
         R"pbdoc(
             Get the current coordinates of the crosslinkers and entanglement links.
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
         &mehp::MEHPForceBalance2::getAtomIdsOfActiveNodes,
         R"pbdoc(
           Get the atom ids of the nodes that are considered active.
           Only crosslink ids are returned (not e.g. entanglement links).

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

                :param direction: The direction in which to compute the active springs (0: x, 1: y, 2: z)
                :param tolerance: springs under this length are considered inactive
           )pbdoc",
         py::arg("direction"),
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_spring",
         &mehp::MEHPForceBalance2::getNrOfActiveSprings,
         R"pbdoc(
           Get the number of active springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_current_spring_vectors",
         &mehp::MEHPForceBalance2::getCurrentSpringDistances,
         R"pbdoc(
           Get the spring vectors.
      )pbdoc")
    .def("get_current_spring_lengths",
         &mehp::MEHPForceBalance2::getCurrentSpringLengths,
         R"pbdoc(
           Get the spring lengths (Euclidean distances of the spring vectors).
      )pbdoc")
    .def("get_overall_strand_lengths",
         &mehp::MEHPForceBalance2::getOverallSpringLengths,
         R"pbdoc(
           Get the sum of the lengths of the springs of each strand.
      )pbdoc")
    .def("get_effective_functionality_of_atoms",
         &mehp::MEHPForceBalance2::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
           Returns the number of active springs connected to each atom, atomId used as index

           :param tolerance: springs under this length are considered inactive
      )pbdoc",
         py::arg("tolerance") = 1e-3)
    .def("get_average_spring_length",
         &mehp::MEHPForceBalance2::getAverageStrandLength,
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
           Returns the universe [of crosslinkers] with the positions of the current state of the simulation.
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
     A quick-and-dirty implementation of the dissipative particle dynamics (DPD) simulation
     with slip-springs as presented by :cite:t:`langeloth_recovering_2013` 
     and :cite:t:`schneider_simulation_2021`.
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
      [](dpd::DPDSimulator& sim,
         const int nSteps,
         const double dt,
         const bool withMC) {
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

          :param file: The file path to the restart file to read
     )pbdoc",
                py::arg("file"))
    .def("config_restart_output",
         &dpd::DPDSimulator::configRestartOutput,
         R"pbdoc(
          Set when to output a restart file.

          Note:
               The filename determines the type of serialization:
               .json, .xml are supported; other file endings will lead to binary serialization (fastest!).

          Caution:
               This method may not be backwards- nor forward-compatible.
               Use the same version of pylimer-tools if you want to be sure that things work.

          :param file: The file path to the restart file to write
          :param output_every: How often to write the restart file (default: 50000)
     )pbdoc",
         py::arg("file"),
         py::arg("output_every") = 50000)
    .def("write_restart_file",
         &dpd::DPDSimulator::writeRestartFile,
         R"pbdoc(
          Explicitly force the writing of a restart file, now!

          :param file: The file path and name of the restart file to be written.
                       Can end in .xml, .json or anything else (-> binary)
     )pbdoc",
         py::arg("file"))
#endif
    .def("config_average_output",
         &dpd::DPDSimulator::configAverageOutput,
         R"pbdoc(
          Set which values to compute averages for.

          :param values: A list of OutputConfiguration structs specifying what to average
     )pbdoc",
         py::arg("values"))
    .def("config_auto_correlator_output",
         &dpd::DPDSimulator::configAutoCorrelatorOutput,
         R"pbdoc(
          Set which values to compute multiple-tau autocorrelation for.
          If you use this, you should cite :cite:t:`ramirez_efficient_2010`.

          :param values: a list of OutputConfiguration structs
          :param num_corr_in: Number of correlations in
          :param p: Parameter p for the autocorrelator
          :param m: Parameter m for the autocorrelator
     )pbdoc",
         py::arg("values"),
         py::arg("num_corr_in") = 32,
         py::arg("p") = 16,
         py::arg("m") = 2)
    .def("config_step_output",
         &dpd::DPDSimulator::configStepOutput,
         R"pbdoc(
          Set which values to log.

          :param values: a list of OutputConfiguration structs
     )pbdoc",
         py::arg("values"))
    .def("config_shift_possibility_empty",
         &dpd::DPDSimulator::configShiftPossibilityEmpty,
         R"pbdoc(
          Configure the possibility of shifting to empty positions.

          This setting affects Monte Carlo moves in the simulation.
         )pbdoc",
         py::arg("shift_possibility_empty") = true)
    .def("config_shift_one_at_a_time",
         &dpd::DPDSimulator::configShiftOneAtATime,
         R"pbdoc(
          Configure whether to shift atoms one at a time.

          This setting affects Monte Carlo move behavior in the simulation.
         )pbdoc",
         py::arg("shift_one_at_a_time") = false)
    .def("config_num_steps_mc",
         &dpd::DPDSimulator::configNumStepsMC,
         R"pbdoc(
          Configure the number of steps to do in one MC sequence.

          :param num_steps: Number of Monte Carlo steps per sequence
     )pbdoc",
         py::arg("num_steps") = 500)
    .def("config_num_steps_dpd",
         &dpd::DPDSimulator::configNumStepsDPD,
         R"pbdoc(
          Configure the number of steps to do in one DPD sequence.

          :param num_steps: Number of DPD steps per sequence
     )pbdoc",
         py::arg("num_steps") = 500)
    .def("config_bond_formation",
         &dpd::DPDSimulator::configBondFormation,
         R"pbdoc(
          Configure how to do bond formation during the run.

          :param num_bonds_to_form (int): The nr of bonds to form in total. Use 0 to stop bond formation.
          :param num_bonds_per_atom_type (dict): The nr of bonds each atom type may have at most (e.g., 2 for strand atoms, 4 for a tertiary crosslinkers)
          :param bond_formation_dist (float): The maximum distance allowed to form bonds
          :param attempt_bond_formation_every (int): attempt to form bonds every this many steps during the simulation run
          :param atom_type_form_from (int): The atom type to start forming bonds from.
          :param atom_type_form_to (int): The atom type to start forming bonds to.
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
          Configure whether a relocation step may happen when a slip-spring has ended at a crosslink.

          Side-effect: if true, the relocations may also happen *to* a slip-spring next to a crosslink.

          :param allow_relocation_in_network (bool): Whether to allow relocation in the network or not.
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

     :param with_slipsprings: Whether to include slip-springs in the returned universe (default: True)
    )pbdoc",
         py::arg("with_slipsprings") = true)
    .def("refresh_current_state",
         &dpd::DPDSimulator::refreshCurrentState,
         R"pbdoc(
          After re-configuring the force-field parameters,
          this method should be called to update the current stress tensor etc.
     )pbdoc")
    .def("get_timestep",
         &dpd::DPDSimulator::getTimestep,
         R"pbdoc(
          Get the timestep used in the simulation.

          :return: The simulation timestep value
         )pbdoc")
    .def("get_current_timestep",
         &dpd::DPDSimulator::getCurrentTimestep,
         R"pbdoc(
          Get the current timestep number.

          :return: The current timestep index
         )pbdoc")
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

          :param file: The file path and name of the restart file to be written.
                       Can end in .xml, .json or anything else (-> binary)
     )pbdoc",
         py::arg("file"))
#endif
    .def("validate_neighbour_list",
         &dpd::DPDSimulator::validateNeighbourlist,
         R"pbdoc(
          Validate the neighbor list consistency.

          :return: True if the neighbor list is valid
         )pbdoc")
    .def("validate_state",
         &dpd::DPDSimulator::validateState,
         R"pbdoc(
          Validate the current simulation state.

          :return: True if the simulation state is valid
         )pbdoc");
}

#endif /* PYBIND_CALC_H */
