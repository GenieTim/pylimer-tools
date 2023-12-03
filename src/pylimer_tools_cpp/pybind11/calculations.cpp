#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/DPDSimulator.h"
#include "../calc/MEHPForceBalance.h"
// #include "../calc/MEHPForceBalance2.h"
#include "../calc/MEHPForceEvaluator.h"
#include "../calc/MEHPForceRelaxation.h"
#include "../calc/MEHPanalysis.h"
#include "../calc/MMTanalysis.h"
#include "../entities/Universe.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;

using namespace pylimer_tools::calc;

namespace pylimer_tools::calc::mehp {
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
    PYBIND11_OVERRIDE_PURE(
      double,                     /* Return type */
      MEHPForceEvaluator,         /* Parent class */
      evaluateStressContribution, /* Name of function in C++ */
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
    const Eigen::VectorXd& u,
    bool requiresGradient) const
  {
    PYBIND11_OVERRIDE_PURE(
      returntype,               /* Return type */
      MEHPForceEvaluator,       /* Parent class */
      evaluateForceSetGradient, /* Name of function in C++ (must match Python
                                   name) */
      n,
      springDistances,
      u,
      requiresGradient /* Argument(s) */
    );
  }

  // actually overriding function, but simplifying for python possibilities
  double evaluateForceSetGradient(const size_t n,
                                  const Eigen::VectorXd& springDistances,
                                  const Eigen::VectorXd& u,
                                  double* grad) const override
  {
    std::pair<double, std::vector<double>> trampolineResult =
      this->evaluateForceAndGradient(n, springDistances, u, grad != nullptr);
    if (grad != nullptr) {
      assert(trampolineResult.second.size() == n);
      for (size_t i = 0; i < n; ++i) {
        grad[i] = trampolineResult.second[i];
      }
    }
    return trampolineResult.first;
  }

  void prepareForEvaluations() override{};
};
}

void
init_pylimer_bound_calc(py::module_& m)
{
  m.def("predictGelationPoint",
        &mmt::predictGelationPoint,
        "Predict the gelation point of a Universe");
  // m.def("computeExtentOfReaction", &mmt::computeExtentOfReaction, "Compute
  // extent of reaction");
  m.def("computeStoichiometricInbalance",
        &mmt::computeStoichiometricInbalance,
        "Compute stoichiometric inbalance");

  /**
   * MEHP
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

  py::class_<mehp::Network>(m, "SimplifiedNetwork", R"pbdoc(
     A more efficient structure of the network for use in MEHP.
     Consists usually only of the cross-links.
 )pbdoc")
    .def_readonly("boxLengths", &mehp::Network::L)
    .def_readonly("volume", &mehp::Network::vol)
    .def_readonly("nrOfNodes", &mehp::Network::nrOfNodes)
    .def_readonly("nrOfSprings", &mehp::Network::nrOfSprings)
    // .def_readonly("nrOfLoops", &mehp::Network::nrOfLoops)
    .def_readonly("coordinates", &mehp::Network::coordinates)
    .def_readonly("oldAtomIds", &mehp::Network::oldAtomIds)
    .def_readonly("springCoordinateIndexA",
                  &mehp::Network::springCoordinateIndexA)
    .def_readonly("springCoordinateIndexB",
                  &mehp::Network::springCoordinateIndexB)
    .def_readonly("springIndexA", &mehp::Network::springIndexA)
    .def_readonly("springIndexB", &mehp::Network::springIndexB)
    // .def_readonly("springIsActive", &mehp::Network::springIsActive)
    ;

  py::class_<mehp::ForceBalanceNetwork>(m, "SimplifiedBalanceNetwork", R"pbdoc(
     A more efficient structure of the network for use in MEHP force balance.
     Consists usually only of the cross- and slip-links.
 )pbdoc")
    .def_readonly("boxLengths", &mehp::ForceBalanceNetwork::L)
    .def_readonly("volume", &mehp::ForceBalanceNetwork::vol)
    .def_readonly("nrOfCrossLinks", &mehp::ForceBalanceNetwork::nrOfNodes)
    .def_readonly("nrOfLinks", &mehp::ForceBalanceNetwork::nrOfLinks)
    .def_readonly("nrOfSprings", &mehp::ForceBalanceNetwork::nrOfSprings)
    .def_readonly("nrOfPartialSprings",
                  &mehp::ForceBalanceNetwork::nrOfPartialSprings)
    // .def_readonly("nrOfLoops", &mehp::Network::nrOfLoops)
    .def_readonly("coordinates", &mehp::ForceBalanceNetwork::coordinates)
    .def_readonly("oldAtomIds", &mehp::ForceBalanceNetwork::oldAtomIds)
    .def_readonly("springCoordinateIndexA",
                  &mehp::ForceBalanceNetwork::springCoordinateIndexA)
    .def_readonly("springCoordinateIndexB",
                  &mehp::ForceBalanceNetwork::springCoordinateIndexB)
    .def_readonly("springPartCoordinateIndexA",
                  &mehp::ForceBalanceNetwork::springPartCoordinateIndexA)
    .def_readonly("springPartCoordinateIndexB",
                  &mehp::ForceBalanceNetwork::springPartCoordinateIndexB)
    .def_readonly("springIndexA", &mehp::ForceBalanceNetwork::springIndexA)
    .def_readonly("springIndexB", &mehp::ForceBalanceNetwork::springIndexB)
    .def_readonly("springPartIndexA",
                  &mehp::ForceBalanceNetwork::springPartIndexA)
    .def_readonly("springPartIndexB",
                  &mehp::ForceBalanceNetwork::springPartIndexB)
    .def_readonly("linkIsSliplink", &mehp::ForceBalanceNetwork::linkIsSliplink)
    .def_readonly("localToGlobalSpringIndex",
                  &mehp::ForceBalanceNetwork::localToGlobalSpringIndex)
    .def_readonly("springIndicesOfLinks",
                  &mehp::ForceBalanceNetwork::springIndicesOfLinks)
    .def_readonly("linkIndicesOfSprings",
                  &mehp::ForceBalanceNetwork::linkIndicesOfSprings)
    .def_readonly("nrOfCrosslinkSwapsEndured",
                  &mehp::ForceBalanceNetwork::nrOfCrosslinkSwapsEndured)
    .def_readonly("springContourLength",
                  &mehp::ForceBalanceNetwork::springsContourLength)
    .def_readonly("oldAtomIds", &mehp::ForceBalanceNetwork::oldAtomIds)
    .def_readonly("partialToFullSpringIndex",
                  &mehp::ForceBalanceNetwork::partialToFullSpringIndex)
    // .def_readonly("springIsActive", &mehp::Network::springIsActive)
    ;

  py::class_<mehp::MEHPForceEvaluator, mehp::PyMEHPForceEvaluator>(
    m, "MEHPForceEvaluator", R"pbdoc(
     The base interface to change the way the force is evaluated during a MEHP run.
    )pbdoc")
    .def(py::init<>())
    .def_property_readonly("network", &mehp::MEHPForceEvaluator::getNetwork)
    .def_property("is2D",
                  &mehp::MEHPForceEvaluator::getIs2D,
                  &mehp::MEHPForceEvaluator::setIs2D)
    //     .def("evaluateForceSetGradient",
    //          py::overload_cast<const size_t,
    //                            const Eigen::VectorXd&,
    //                            const Eigen::VectorXd&,
    //                            bool>(
    //            &mehp::MEHPForceEvaluator::evaluateForceSetGradient))
    .def("evaluateStressContribution",
         &mehp::MEHPForceEvaluator::evaluateStressContribution,
         R"pbdoc(
          An evaluation of the stress-contribution.

          :param springDistances: the three coordinate differences for one spring.
          :param i: the row index of the stress tensor
          :param j: the column index of the stress tensor
    )pbdoc",
         py::arg("springDistances"),
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
  //      :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceEvaluator.evaluateStressContribution`.

  //      :param n: the dimensionality of the problem (the nr. of spring
  //      coordinates) :param springDistances: the sequential (x, y, z) spring
  //      distances :param displacements: the displacements from the original
  //      coordinates
  //           (accessible by
  //           :func:`~pylimer_tools_cpp.pylimer_tools_cpp.CustomMEHPForceEvaluator.getNetwork().coordinates`)
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
    m, "SimpleSpringMEHPForceEvaluator", R"pbdoc(
     This is equal to a spring evaluator for Gaussian chains.

     The force for a certain spring is given by:
     :math:`f = 0.5 \cdot \kappa r`, 
     where :math:`r` is the spring [between cross-linkers] length.

     Recommended optimization algorithm: "LD_LBFGS"

     :param kappa: the spring constant :math:`\kappa`
    )pbdoc")
    .def(py::init<double>(), py::arg("kappa") = 1.0);

  py::class_<mehp::NonGaussianSpringForceEvaluator, mehp::MEHPForceEvaluator>(
    m, "NonGaussianSpringForceEvaluator", R"pbdoc(
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

  py::class_<mehp::MEHPForceRelaxation>(m, "MEHPForceRelaxation", R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
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
          :param crosslinkerType: The atom type of the cross-linkers. Needed to reduce the network.
          :param is2D: Whether to ignore the z direction.
          :param forceEvaluator: The force evaluator to use
          :param kappa: The spring constant
          :param remove2functionalCrosslinkers: Whether to replace two-functional cross-links with a "normal" chain bead
          :param removeDanglingChains: Whether to remove dangling chains before running the simulation. 
               **Caution*: Removing the dangling chains will result in incorrect results fo the computation of 
               :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getSolubleWeightFraction()` and
               :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getDanglingWeightFraction()`
          )pbdoc",
         py::arg("universe"),
         py::arg("crosslinkerType") = 2,
         py::arg("is2D") = false,
         py::arg("forceEvaluator") = nullptr,
         py::arg("kappa") = 1.0,
         py::arg("remove2functionalCrosslinkers") = true,
         py::arg("removeDanglingChains") = false)
    .def("runForceRelaxation",
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
         py::arg("maxNrOfSteps") = 250000,
         py::arg("xTolerance") = 1e-12,
         py::arg("fTolerance") = 1e-9)
    // .def("getForceEvaluator", &mehp::MEHPForceRelaxation::getForceEvaluator,
    // R"pbdoc(
    //      Query the currently used force evaluator.
    // )pbdoc")
    .def("setForceEvaluator",
         &mehp::MEHPForceRelaxation::setForceEvaluator,
         R"pbdoc(
          Reset the currently used force evaluator.
     )pbdoc")
    .def("getForce",
         &mehp::MEHPForceRelaxation::getForce,
         R"pbdoc(
          Returns the force at the current state of the simulation.
     )pbdoc")
    .def("getResiduals",
         &mehp::MEHPForceRelaxation::getResiduals,
         R"pbdoc(
          Returns the residuals at the current state of the simulation.
     )pbdoc")
    .def("getResidualNorm",
         &mehp::MEHPForceRelaxation::getResidualNorm,
         R"pbdoc(
          Returns the residual norm at the current state of the simulation.
     )pbdoc")
    .def("getPressure",
         &mehp::MEHPForceRelaxation::getPressure,
         R"pbdoc(
          Returns the pressure at the current state of the simulation.
     )pbdoc")
    .def("getStressTensor",
         &mehp::MEHPForceRelaxation::getStressTensor,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc")
    .def("getGammaFactor",
         &mehp::MEHPForceRelaxation::getGammaFactor,
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
          
          :param r0squared: The denominator in the equation of :math:`\Gamma`. If :math:`-1.0` (default), the network is used for determination (which is not accurate). For phantom systems, the correct value is :math:`Nb^2`.
               For other systems, the value could be determined by `~pylimer_tools_cpp.pylimer_tools_cpp.Universe.computeMeanEndToEndDistance` on the melt system.
          :param nrOfChains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains. 
     )pbdoc",
         py::arg("r0squared") = -1.0,
         py::arg("nrOfChains") = -1)
    .def("getNrOfNodes", &mehp::MEHPForceRelaxation::getNrOfNodes, R"pbdoc(
           Get the number of nodes considered in this simulation.
     )pbdoc")
    .def("getNrOfSprings",
         &mehp::MEHPForceRelaxation::getNrOfSprings,
         R"pbdoc(
          Get the number of springs considered in this simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc")
    .def("getIdsOfActiveNodes",
         &mehp::MEHPForceRelaxation::getIdsOfActiveNodes,
         R"pbdoc(
          Get the atom ids of the nodes that are considered active.

          :param tolerance: springs under this length are considered inactive. A node is active if it has > 2 active springs.
          :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
          :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
               Use a value < 0 to indicate that there is no maximum number of active connections.
     )pbdoc",
         py::arg("tolerance") = 0.1,
         py::arg("minimumNrOfActiveConnections") = 2,
         py::arg("maximumNrOfActiveConnections") = -1)
    .def("getNrOfActiveNodes",
         &mehp::MEHPForceRelaxation::getNrOfActiveNodes,
         R"pbdoc(
           Get the number of active nodes remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive.
          :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
          :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
               Use a value < 0 to indicate that there is no maximum number of active connections.
     )pbdoc",
         py::arg("tolerance") = 0.1,
         py::arg("minimumNrOfActiveConnections") = 2,
         py::arg("maximumNrOfActiveConnections") = -1)
    .def("getNrOfActiveSprings",
         &mehp::MEHPForceRelaxation::getNrOfActiveSprings,
         R"pbdoc(
           Get the number of active springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getSolubleWeightFraction",
         &mehp::MEHPForceRelaxation::getSolubleWeightFraction,
         R"pbdoc(
          Compute the weight fraction of springs connected to active
          springs (any depth). 
          
          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getDanglingWeightFraction",
         &mehp::MEHPForceRelaxation::getDanglingWeightFraction,
         R"pbdoc(
          Compute the weight fraction of non-active springs

          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getEffectiveFunctionalityOfAtoms",
         &mehp::MEHPForceRelaxation::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
          Returns the number of active springs connected to each atom, atomId used as index

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def(
      "getSpringLengths", &mehp::MEHPForceRelaxation::getSpringLengths, R"pbdoc(
          Get the current lengths for all the springs.

          Returns:
               - distances: a vector of size nrOfSprings, with each the norm of the distances
     )pbdoc")
    .def("getSpringDistances",
         &mehp::MEHPForceRelaxation::getSpringDistances,
         R"pbdoc(
          Get the current coordinate differences for all the springs.

          Returns:
               - distances: a vector of size 3*nrOfSprings, with each x, y, z values of the springs
     )pbdoc")
    .def("getAverageSpringLength",
         &mehp::MEHPForceRelaxation::getAverageSpringLength,
         R"pbdoc(
           Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()`,
           this value is normalized by the number of springs rather than the number of chains.
     )pbdoc")
    .def("getDefaultR0Square",
         &mehp::MEHPForceRelaxation::getDefaultR0Square,
         R"pbdoc(
           Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()` for :math:`\langle R_{0,\eta}^2\rangle`.
     )pbdoc")
    .def("getDefaultNrOfChains",
         &mehp::MEHPForceRelaxation::getDefaultNrOfChains,
         R"pbdoc(
          Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()` for normalizing the distances.`.
     )pbdoc")
    .def("getCurrentDisplacements",
         &mehp::MEHPForceRelaxation::getCurrentDisplacements,
         R"pbdoc(
          Returns the current displacement.
     )pbdoc")
    .def("getNrOfIterations",
         &mehp::MEHPForceRelaxation::getNrOfIterations,
         R"pbdoc(
          Returns the number of iterations used for force relaxation.
     )pbdoc")
    .def("getExitReason", &mehp::MEHPForceRelaxation::getExitReason, R"pbdoc(
           Returns the reason for termination of the simulation
     )pbdoc")
    .def("getCrosslinkerVerse",
         &mehp::MEHPForceRelaxation::getCrosslinkerVerse,
         R"pbdoc(
          Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
     )pbdoc");

  py::enum_<mehp::BalanceRunMode>(m, "BalanceRunMode")
    .value("EIGEN_ALL", mehp::BalanceRunMode::EIGEN_ALL)
    .value("EIGEN_RANDOM", mehp::BalanceRunMode::EIGEN_RANDOM)
    .value("EIGEN_HEURISTIC", mehp::BalanceRunMode::EIGEN_HEURISTIC)
    .value("EIGEN_STRANDS", mehp::BalanceRunMode::EIGEN_STRANDS)
    .value("ITERATIVE", mehp::BalanceRunMode::ITERATIVE);

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
    .value("ALL_MC_CYCLE", mehp::LinkSwappingMode::ALL_MC_CYCLE);

  py::class_<mehp::MEHPForceBalance>(m, "MEHPForceBalance", R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
     )pbdoc")
    .def(py::init<pe::Universe, int, bool, double, bool>(),
         R"pbdoc(
          Instantiate the simulator for a certain universe.

          :param universe: the universe to simulate with
          :param crosslinkerType: The atom type of the cross-linkers. Needed to reduce the network.
          :param is2D: Whether to ignore the z direction.
          :param forceEvaluator: The force evaluator to use
          )pbdoc",
         py::arg("universe"),
         py::arg("crosslinkerType") = 2,
         py::arg("is2D") = false,
         py::arg("kappa") = 1.0,
         py::arg("remove2functionalCrosslinkers") = true)
    .def("__copy__",
         [](const mehp::MEHPForceBalance& self) {
           return mehp::MEHPForceBalance(self);
         })
    .def_property_readonly("network", &mehp::MEHPForceBalance::getNetwork)
    //     .def("validateNetwork",
    //          py::overload_cast<>(&mehp::MEHPForceBalance::validateNetwork),
    //          R"pbdoc(
    //            Validates the internal structures.
    //      )pbdoc")
    .def("runForceRelaxation",
         &mehp::MEHPForceBalance::runForceRelaxation,
         R"pbdoc(
          Run the simulation.
          Note that the final state of the minimization is persisted and reused if you use this method again.
          This is useful if you want to run a global optimization first and add a local one afterwards.
          As a consequence though, you cannot simply benchmark only this method; you must include the setup.

          :param runMode: Choice of the mode to run the simulation with.
          :param damping: For certain run modes, a damping factor helps to improve performance.
          :param maxNrOfSteps: The maximum number of steps to do during the simulation.
          :param xTolerance: The tolerance of the displacements as an exit condition.
          :param innerMaxNrOfSteps: The maximum number of steps to do per iteration during the slip-link displacements.
          :param innerXTolerance: The tolerance of the displacements of the slip-link as an inner exit condition.
          :param innerAlphaTolerance: The tolerance of the contour-length when slipping the slip-link as an inner exit condition.
          )pbdoc",
         py::arg("runMode") = mehp::BalanceRunMode::ITERATIVE,
         py::arg("damping") = 1.0,
         py::arg("maxNrOfSteps") = 250000,
         py::arg("xTolerance") = 1e-12,
         py::arg("initialResidualNorm") = -1.0,
         py::arg("simplificationMode") =
           mehp::StructureSimplificationMode::NO_SIMPLIFICATION,
         py::arg("inactiveRemovalCutoff") = -1.0,
         py::arg("outputFrequency") = 50,
         py::arg("doInnerIterations") = false,
         py::arg("allowSlipLinksToPassEachOther") =
           mehp::LinkSwappingMode::NO_SWAPPING,
         py::arg("swappingFrequency") = 10,
         py::arg("oneOverSpringPartitionUpperLimit") = 1.0,
         py::arg("nrOfCrosslinkSwapsAllowedPerSliplink") = -1)
    .def("deformTo",
         &mehp::MEHPForceBalance::deformTo,
         R"pbdoc()pbdoc",
         py::arg("newBox"))
    .def("inspectLinkDisplacementToMeanPositionUpdate",
         &mehp::MEHPForceBalance::inspectLinkDisplacementToMeanPositionUpdate,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("damping") = 1.0)
    .def("swapSlipLinksInclXlinks",
         &mehp::MEHPForceBalance::swapSlipLinksInclXlinks)
    .def("moveSlipLinksToTheirBestBranch",
         &mehp::MEHPForceBalance::moveSlipLinksToTheirBestBranch)
    .def("getForceOn",
         &mehp::MEHPForceBalance::getForceOn,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
    .def("inspectDisplacementToMeanPositionUpdate",
         &mehp::MEHPForceBalance::inspectDisplacementToMeanPositionUpdate,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
    .def("inspectSpringPartitionUpdate",
         &mehp::MEHPForceBalance::inspectSpringPartitionUpdate,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"))
    .def("inspectParametrisationOptimsationForLink",
         &mehp::MEHPForceBalance::inspectParametrisationOptimsationForLink,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("displacements"),
         py::arg("springPartitions"),
         py::arg("maxNrOfSteps") = 100,
         py::arg("alpha_tol") = 1e-9,
         py::arg("minNrOfSteps") = 1,
         py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
    .def("getSpringpartitionIndicesOfSliplink",
         &mehp::MEHPForceBalance::getSpringpartitionIndicesOfSliplink,
         R"pbdoc()pbdoc",
         py::arg("network"),
         py::arg("linkIdx"))
    .def("getNeighbourLinkIndices",
         &mehp::MEHPForceBalance::getNeighbourLinkIndices,
         R"pbdoc()pbdoc",
         py::arg("network"),
         py::arg("linkIdx"))
    .def("evaluateDistanceBetween",
         &mehp::MEHPForceBalance::evaluateDistanceBetween,
         R"pbdoc()pbdoc",
         py::arg("network"),
         py::arg("displacements"),
         py::arg("linkIndexA"),
         py::arg("linkIndexB"),
         py::arg("is2D") = false)
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
    .def("getPressure",
         &mehp::MEHPForceBalance::getPressure,
         R"pbdoc(
          Returns the pressure at the current state of the simulation.
     )pbdoc")
    .def("getSolubleWeightFraction",
         &mehp::MEHPForceBalance::getSolubleWeightFraction,
         R"pbdoc(
          Compute the weight fraction of springs connected to active
          springs (any depth). 
          
          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getDanglingWeightFraction",
         &mehp::MEHPForceBalance::getDanglingWeightFraction,
         R"pbdoc(
          Compute the weight fraction of non-active springs

          Caution: ignores atom masses.
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("addSlipLinks",
         py::overload_cast<const std::vector<size_t>&,
                           const std::vector<size_t>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const bool>(&mehp::MEHPForceBalance::addSlipLinks),
         R"pbdoc(
          Add the slip-links
     )pbdoc",
         py::arg("strandIdx1"),
         py::arg("strandIdx2"),
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("alpha1"),
         py::arg("alpha2"),
         py::arg("clampAlpha") = false)
    .def("randomlyAddSlipLinks",
         &mehp::MEHPForceBalance::randomlyAddSliplinks,
         R"pbdoc()pbdoc",
         py::arg("nrOfSlipLinksToSample"),
         py::arg("cutoff") = 2.0,
         py::arg("minimumNrOfSliplinks") = 0,
         py::arg("sameStrandCutoff") = 2,
         py::arg("excludeCrosslinks") = false,
         py::arg("seed") = -1)
    .def("addSliplinksBasedOnCycles",
         &mehp::MEHPForceBalance::addSliplinksBasedOnCycles,
         R"pbdoc()pbdoc",
         py::arg("maxLoopLength") = -1)
    .def("getStressTensor",
         &mehp::MEHPForceBalance::getStressTensorArray,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc",
         py::arg("oneOverSpringPartitionUpperLimit") = 1.)
    .def("getStressTensorLinkBased",
         &mehp::MEHPForceBalance::getStressTensorArrayLinkBased,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc",
         py::arg("oneOverSpringPartitionUpperLimit") = 1.,
         py::arg("xlinksOnly") = false)
    .def("getGammaFactor",
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
          
          :param r0squared: The denominator in the equation of :math:`\Gamma`. If :math:`-1.0` (default), the network is used for determination (which is not accurate). For phantom systems, the correct value is :math:`Nb^2`.
               For other systems, the value could be determined by `~pylimer_tools_cpp.pylimer_tools_cpp.Universe.computeMeanEndToEndDistance` on the melt system.
          :param nrOfChains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains. 
     )pbdoc",
         py::arg("r0squared") = -1.0,
         py::arg("nrOfChains") = -1)
    .def("getNrOfNodes", &mehp::MEHPForceBalance::getNrOfNodes, R"pbdoc(
           Get the number of nodes considered in this simulation.
     )pbdoc")
    .def("getNrOfSprings",
         &mehp::MEHPForceBalance::getNrOfSprings,
         R"pbdoc(
          Get the number of springs considered in this simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc")
    .def("getSpringPartitions",
         &mehp::MEHPForceBalance::getSpringPartitions,
         R"pbdoc(
          Get the current spring partitions.
     )pbdoc")
    .def("setSpringPartitions",
         &mehp::MEHPForceBalance::setSpringPartitions,
         R"pbdoc(
          Set the current spring partitions.
     )pbdoc")
    .def("getDisplacements",
         &mehp::MEHPForceBalance::getCurrentDisplacements,
         R"pbdoc(
          Get the current link displacements.
     )pbdoc")
    .def("setDisplacements",
         &mehp::MEHPForceBalance::setCurrentDisplacements,
         R"pbdoc(
          Set the current link displacements.
     )pbdoc")
    .def("setSpringContourLengths",
         &mehp::MEHPForceBalance::setSpringContourLengths,
         R"pbdoc(
          Set/overwrite the contour lengths.
     )pbdoc")
    .def("getDisplacementResidualNorm",
         &mehp::MEHPForceBalance::getDisplacementResidualNorm,
         R"pbdoc(
          Get the current link displacement residual norm.
     )pbdoc")
    .def("getIdsOfActiveNodes",
         &mehp::MEHPForceBalance::getIdsOfActiveNodes,
         R"pbdoc(
          Get the atom ids of the nodes that are considered active.

          :param tolerance: springs under this length are considered inactive. A node is active if it has > 2 active springs.
          :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
          :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
               Use a value < 0 to indicate that there is no maximum number of active connections.
     )pbdoc",
         py::arg("tolerance") = 0.1,
         py::arg("minimumNrOfActiveConnections") = 2,
         py::arg("maximumNrOfActiveConnections") = -1,
         py::arg("usePartial") = false)
    .def("getNrOfActiveNodes",
         &mehp::MEHPForceBalance::getNrOfActiveNodes,
         R"pbdoc(
           Get the number of active nodes remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive.
          :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
          :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
               Use a value < 0 to indicate that there is no maximum number of active connections.
          :param usePartial: Whether to use the partial spring distances rather than the total (set to true if you want primary loop contributors)
     )pbdoc",
         py::arg("tolerance") = 0.1,
         py::arg("minimumNrOfActiveConnections") = 2,
         py::arg("maximumNrOfActiveConnections") = -1,
         py::arg("usePartial") = false)
    .def("getNrOfActiveSprings",
         &mehp::MEHPForceBalance::getNrOfActiveSprings,
         R"pbdoc(
           Get the number of active springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getNrOfActivePartialSprings",
         &mehp::MEHPForceBalance::getNrOfActivePartialSprings,
         R"pbdoc(
           Get the number of active partial springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getEffectiveFunctionalityOfAtoms",
         &mehp::MEHPForceBalance::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
          Returns the number of active springs connected to each atom, atomId used as index

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getAverageSpringLength",
         &mehp::MEHPForceBalance::getAverageSpringLength,
         R"pbdoc(
           Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()`,
           this value is normalized by the number of springs rather than the number of chains.
     )pbdoc")
    .def("getDefaultR0Square",
         &mehp::MEHPForceBalance::getDefaultR0Square,
         R"pbdoc(
           Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()` for :math:`\langle R_{0,\eta}^2\rangle`.
     )pbdoc")
    .def("getDefaultNrOfChains",
         &mehp::MEHPForceBalance::getDefaultNrOfChains,
         R"pbdoc(
          Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()` for normalizing the distances.`.
     )pbdoc")
    .def("getNrOfIterations",
         &mehp::MEHPForceBalance::getNrOfIterations,
         R"pbdoc(
          Returns the number of iterations used for force relaxation.
     )pbdoc")
    .def("getExitReason", &mehp::MEHPForceBalance::getExitReason, R"pbdoc(
           Returns the reason for termination of the simulation
     )pbdoc")
    .def("getCrosslinkerVerse",
         &mehp::MEHPForceBalance::getCrosslinkerVerse,
         R"pbdoc(
          Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
     )pbdoc");

  ////////////////////////////////////////////////////////////////

  //   py::class_<mehp::MEHPForceBalance2>(m, "MEHPForceBalance2", R"pbdoc(
  //     A small simulation tool for quickly minimizing the force between the
  //     cross-linker beads.
  //      )pbdoc")
  //     .def(py::init<pe::Universe, int, bool>(),
  //          R"pbdoc(
  //           Instantiate the simulator for a certain universe.

  //           :param universe: the universe to simulate with
  //           :param crosslinkerType: The atom type of the cross-linkers.
  //           Needed to reduce the network. :param is2D: Whether to ignore the
  //           z direction. :param forceEvaluator: The force evaluator to use
  //           )pbdoc",
  //          py::arg("universe"),
  //          py::arg("crosslinkerType") = 2,
  //          py::arg("is2D") = false)
  //     .def_property_readonly("network", &mehp::MEHPForceBalance2::getNetwork)
  //     .def("validateNetwork",
  //          py::overload_cast<>(&mehp::MEHPForceBalance2::validateNetwork),
  //          R"pbdoc(
  //            Validates the internal structures.
  //      )pbdoc")
  //     .def("runForceRelaxation",
  //          &mehp::MEHPForceBalance2::runForceRelaxation,
  //          R"pbdoc(
  //           Run the simulation.
  //           Note that the final state of the minimization is persisted and
  //           reused if you use this method again. This is useful if you want
  //           to run a global optimization first and add a local one
  //           afterwards. As a consequence though, you cannot simply benchmark
  //           only this method; you must include the setup. )pbdoc",
  //          py::arg("algorithm") = "LD_MMA",
  //          py::arg("maxNrOfSteps") = 250000,
  //          py::arg("xTolerance") = 1e-12,
  //          py::arg("fTolerance") = 1e-9,
  //          py::arg("constraintTol") = 1e-9)
  //     .def("inspectLinkDisplacementToMeanPositionUpdate",
  //          &mehp::MEHPForceBalance2::inspectLinkDisplacementToMeanPositionUpdate,
  //          R"pbdoc()pbdoc",
  //          py::arg("linkIdx"),
  //          py::arg("damping") = 1.0)
  //     .def("getForceOn",
  //          &mehp::MEHPForceBalance2::getForceOn,
  //          R"pbdoc()pbdoc",
  //          py::arg("linkIdx"),
  //          py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
  //     .def("inspectDisplacementToMeanPositionUpdate",
  //          &mehp::MEHPForceBalance2::inspectDisplacementToMeanPositionUpdate,
  //          R"pbdoc()pbdoc",
  //          py::arg("linkIdx"),
  //          py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
  //     .def("inspectSpringPartitionUpdate",
  //          &mehp::MEHPForceBalance2::inspectSpringPartitionUpdate,
  //          R"pbdoc()pbdoc",
  //          py::arg("linkIdx"))
  //     .def("inspectParametrisationOptimsationForLink",
  //          &mehp::MEHPForceBalance2::inspectParametrisationOptimsationForLink,
  //          R"pbdoc()pbdoc",
  //          py::arg("linkIdx"),
  //          py::arg("displacements"),
  //          py::arg("springPartitions"),
  //          py::arg("maxNrOfSteps") = 100,
  //          py::arg("alpha_tol") = 1e-9,
  //          py::arg("distanceBackTolerance") = 1e-9,
  //          py::arg("residualNormSTolerance") = 1e-20,
  //          py::arg("minNrOfSteps") = 1,
  //          py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
  //     .def("getSpringpartitionIndicesOfSliplink",
  //          &mehp::MEHPForceBalance2::getSpringpartitionIndicesOfSliplink,
  //          R"pbdoc()pbdoc",
  //          py::arg("linkIdx"))
  //     // .def("getForceEvaluator",
  //     &mehp::MEHPForceBalance2::getForceEvaluator,
  //     // R"pbdoc(
  //     //      Query the currently used force evaluator.
  //     // )pbdoc")
  //     //     .def("setForceEvaluator",
  //     //          &mehp::MEHPForceBalance2::setForceEvaluator,
  //     //          R"pbdoc(
  //     //           Reset the currently used force evaluator.
  //     //      )pbdoc")
  //     //     .def("getForce",
  //     //          &mehp::MEHPForceBalance2::getForce,
  //     //          R"pbdoc(
  //     //           Returns the force at the current state of the simulation.
  //     //      )pbdoc")
  //     //     .def("getResidualNorm",
  //     //          &mehp::MEHPForceBalance2::getResidualNorm,
  //     //          R"pbdoc(
  //     //           Returns the residual norm at the current state of the
  //     //           simulation.
  //     //      )pbdoc")
  //     .def("getPressure",
  //          &mehp::MEHPForceBalance2::getPressure,
  //          R"pbdoc(
  //           Returns the pressure at the current state of the simulation.
  //      )pbdoc")
  //     .def("addSlipLinks",
  //          py::overload_cast<const std::vector<size_t>,
  //                            const std::vector<size_t>,
  //                            const std::vector<double>,
  //                            const std::vector<double>,
  //                            const std::vector<double>,
  //                            const std::vector<double>,
  //                            const std::vector<double>,
  //                            const
  //                            bool>(&mehp::MEHPForceBalance2::addSlipLinks),
  //          R"pbdoc(
  //           Add the slip-links
  //      )pbdoc",
  //          py::arg("strandIdx1"),
  //          py::arg("strandIdx2"),
  //          py::arg("x"),
  //          py::arg("y"),
  //          py::arg("z"),
  //          py::arg("alpha1"),
  //          py::arg("alpha2"),
  //          py::arg("clampAlpha") = false)
  //     .def("getStressTensor",
  //          &mehp::MEHPForceBalance2::getStressTensor,
  //          R"pbdoc(
  //           Returns the stress tensor at the current state of the simulation.
  //      )pbdoc",
  //          py::arg("oneOverSpringPartitionUpperLimit") = -1)
  //     .def("getStressTensorLinkBased",
  //          &mehp::MEHPForceBalance2::getStressTensorLinkBased,
  //          R"pbdoc(
  //           Returns the stress tensor at the current state of the simulation.
  //      )pbdoc",
  //          py::arg("oneOverSpringPartitionUpperLimit") = -1,
  //          py::arg("xlinksOnly") = false)
  //     .def("getGammaFactor",
  //          &mehp::MEHPForceBalance2::getGammaFactor,
  //          R"pbdoc(
  //           Computes the gamma factor as part of the ANT/MEHP formulism,
  //           i.e.:

  //           :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with
  //           :math:`\gamma_{\eta} = \frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
  //           which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
  //           where :math:`\eta` is the index of a particular strand,
  //           :math:`R_{0}^2` is the melt mean square end to end distance, in
  //           phantom systems :math:`$= N_{\eta}*b^2$` :math:`N_{\eta}` is the
  //           number of atoms in this strand :math:`\eta`, :math:`b` its mean
  //           square bond length, :math:`T` the temperature and :math:`k_B`
  //           Boltzmann's constant.

  //           :param r0squared: The denominator in the equation of
  //           :math:`\Gamma`. If :math:`-1.0` (default), the network is used
  //           for determination (which is not accurate). For phantom systems,
  //           the correct value is :math:`Nb^2`.
  //                For other systems, the value could be determined by
  //                `~pylimer_tools_cpp.pylimer_tools_cpp.Universe.computeMeanEndToEndDistance`
  //                on the melt system.
  //           :param nrOfChains: the value to normalize the sum of square
  //           distances by. Usually (and default if :math:`< 0`) the nr of
  //           chains.
  //      )pbdoc",
  //          py::arg("r0squared") = -1.0,
  //          py::arg("nrOfChains") = -1)
  //     .def("getNrOfNodes", &mehp::MEHPForceBalance2::getNrOfNodes, R"pbdoc(
  //            Get the number of nodes considered in this simulation.
  //      )pbdoc")
  //     .def("getNrOfSprings",
  //          &mehp::MEHPForceBalance2::getNrOfSprings,
  //          R"pbdoc(
  //           Get the number of springs considered in this simulation.

  //           :param tolerance: springs under this length are considered
  //           inactive
  //      )pbdoc")
  //     .def("getSpringPartitions",
  //          &mehp::MEHPForceBalance2::getSpringPartitions,
  //          R"pbdoc(
  //           Get the current spring partitions.
  //      )pbdoc")
  //     .def("setSpringPartitions",
  //          &mehp::MEHPForceBalance2::setSpringPartitions,
  //          R"pbdoc(
  //           Set the current spring partitions.
  //      )pbdoc")
  //     .def("getDisplacements",
  //          &mehp::MEHPForceBalance2::getCurrentDisplacements,
  //          R"pbdoc(
  //           Get the current link displacements.
  //      )pbdoc")
  //     .def("setDisplacements",
  //          &mehp::MEHPForceBalance2::setCurrentDisplacements,
  //          R"pbdoc(
  //           Set the current link displacements.
  //      )pbdoc")
  //     .def("setSpringContourLengths",
  //          &mehp::MEHPForceBalance2::setSpringContourLengths,
  //          R"pbdoc(
  //           Set/overwrite the contour lengths.
  //      )pbdoc")
  //     .def("getDisplacementResidualNorm",
  //          py::overload_cast<double>(
  //            &mehp::MEHPForceBalance2::getDisplacementResidualNorm),
  //          R"pbdoc(
  //           Get the current link displacement residual norm.
  //      )pbdoc")
  //     .def("getIdsOfActiveNodes",
  //          &mehp::MEHPForceBalance2::getIdsOfActiveNodes,
  //          R"pbdoc(
  //           Get the atom ids of the nodes that are considered active.

  //           :param tolerance: springs under this length are considered
  //           inactive. A node is active if it has > 2 active springs. :param
  //           minimumNrOfActiveConnections:  A node is active if it has equal
  //           or more than this number of active springs. :param
  //           maximumNrOfActiveConnections:  A node is active if it has equal
  //           or less than this number of active springs.
  //                Use a value < 0 to indicate that there is no maximum number
  //                of active connections.
  //      )pbdoc",
  //          py::arg("tolerance") = 0.1,
  //          py::arg("minimumNrOfActiveConnections") = 2,
  //          py::arg("maximumNrOfActiveConnections") = -1,
  //          py::arg("usePartial") = false)
  //     .def("getNrOfActiveNodes",
  //          &mehp::MEHPForceBalance2::getNrOfActiveNodes,
  //          R"pbdoc(
  //            Get the number of active nodes remaining after running the
  //            simulation.

  //           :param tolerance: springs under this length are considered
  //           inactive. :param minimumNrOfActiveConnections:  A node is active
  //           if it has equal or more than this number of active springs.
  //           :param maximumNrOfActiveConnections:  A node is active if it has
  //           equal or less than this number of active springs.
  //                Use a value < 0 to indicate that there is no maximum number
  //                of active connections.
  //           :param usePartial: Whether to use the partial spring distances
  //           rather than the total (set to true if you want primary loop
  //           contributors)
  //      )pbdoc",
  //          py::arg("tolerance") = 0.1,
  //          py::arg("minimumNrOfActiveConnections") = 2,
  //          py::arg("maximumNrOfActiveConnections") = -1,
  //          py::arg("usePartial") = false)
  //     .def("getNrOfActiveSprings",
  //          &mehp::MEHPForceBalance2::getNrOfActiveSprings,
  //          R"pbdoc(
  //            Get the number of active springs remaining after running the
  //            simulation.

  //           :param tolerance: springs under this length are considered
  //           inactive
  //      )pbdoc",
  //          py::arg("tolerance") = 0.1)
  //     .def("getNrOfActivePartialSprings",
  //          &mehp::MEHPForceBalance2::getNrOfActivePartialSprings,
  //          R"pbdoc(
  //            Get the number of active partial springs remaining after running
  //            the simulation.

  //           :param tolerance: springs under this length are considered
  //           inactive
  //      )pbdoc",
  //          py::arg("tolerance") = 0.1)
  //     .def("getEffectiveFunctionalityOfAtoms",
  //          &mehp::MEHPForceBalance2::getEffectiveFunctionalityOfAtoms,
  //          R"pbdoc(
  //           Returns the number of active springs connected to each atom,
  //           atomId used as index

  //           :param tolerance: springs under this length are considered
  //           inactive
  //      )pbdoc",
  //          py::arg("tolerance") = 0.1)
  //     .def("getAverageSpringLength",
  //          &mehp::MEHPForceBalance2::getAverageSpringLength,
  //          R"pbdoc(
  //            Get the average length of the springs. Note that in contrast to
  //            :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()`,
  //            this value is normalized by the number of springs rather than
  //            the number of chains.
  //      )pbdoc")
  //     .def("getDefaultR0Square",
  //          &mehp::MEHPForceBalance2::getDefaultR0Square,
  //          R"pbdoc(
  //            Returns the value effectively used in
  //            :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()`
  //            for :math:`\langle R_{0,\eta}^2\rangle`.
  //      )pbdoc")
  //     .def("getDefaultNrOfChains",
  //          &mehp::MEHPForceBalance2::getDefaultNrOfChains,
  //          R"pbdoc(
  //           Returns the value effectively used in
  //           :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()`
  //           for normalizing the distances.`.
  //      )pbdoc")
  //     .def("getNrOfIterations",
  //          &mehp::MEHPForceBalance2::getNrOfIterations,
  //          R"pbdoc(
  //           Returns the number of iterations used for force relaxation.
  //      )pbdoc")
  //     .def("getExitReason", &mehp::MEHPForceBalance2::getExitReason, R"pbdoc(
  //            Returns the reason for termination of the simulation
  //      )pbdoc")
  //     .def("getCrosslinkerVerse",
  //          &mehp::MEHPForceBalance2::getCrosslinkerVerse,
  //          R"pbdoc(
  //           Returns the universe [of cross-linkers] with the positions of the
  //           current state of the simulation.
  //      )pbdoc",
  //          py::arg("newCrosslinkerType") = 2);

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
    .value("MEAN_B", ComputedDoubleValues::MEAN_B)
    .value("MAX_B", ComputedDoubleValues::MAX_B)
    .value("MSD", ComputedDoubleValues::MSD);

  py::enum_<ComputedIntValues>(m, "ComputedIntValues")
    .value("STEP", ComputedIntValues::STEP)
    .value("NUM_SHIFT", ComputedIntValues::NUM_SHIFT)
    .value("NUM_RELOC", ComputedIntValues::NUM_RELOC);

  py::class_<OutputConfiguration>(m, "OutputConfiguration")
    .def(py::init<>(), "Get an instance of this struct")
    .def_readwrite("intValues", &OutputConfiguration::intValues)
    .def_readwrite("doubleValues", &OutputConfiguration::doubleValues)
    .def_readwrite(
      "filename",
      &OutputConfiguration::filename,
      R"pbdoc(The file to write to. Empty means standard output (console).)pbdoc")
    .def_readwrite("outputEvery",
                   &OutputConfiguration::outputEvery,
                   R"pbdoc(How often to write the values to the output. 
      For averages, this value also says how many values will be averaged.
     )pbdoc");

  /**
   * DPD Simulations
   */
  py::class_<dpd::DPDSimulator>(m,
                                "DPDSimulator",
                                R"pbdoc(
          A quick-and-dirty implementation of the DPD simulation
          with slip-springs as presented by Langeloth et al.
     )pbdoc")
    .def(
      py::init<const pe::Universe, const int, const bool, const std::string>(),
      "Get an instance of this class",
      py::arg("universe"),
      py::arg("crosslinker_type") = 2,
      py::arg("is_2D") = false,
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
      "runSimulation",
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
    .def("createSlipSprings",
         &dpd::DPDSimulator::createSlipSprings,
         R"pbdoc(
          Randomly add the specified number of slip-springs to neighbours within the specified cut-offs.
     )pbdoc",
         py::arg("num"),
         py::arg("bond_type") = 9)
    .def("configA",
         &dpd::DPDSimulator::configA,
         R"pbdoc(
          Configure the force-field (pair-style) parameter `A`.
     )pbdoc",
         py::arg("A") = 25.)
    .def("configSigma",
         &dpd::DPDSimulator::configSigma,
         R"pbdoc(
          Configure the force-field (pair-style) parameter `\sigma`.
     )pbdoc",
         py::arg("sigma") = 3.)
    .def("configSpringConstant",
         &dpd::DPDSimulator::configSpringConstant,
         R"pbdoc(
          Configure the force-field (bond-style) parameter `k`, the spring constant.
     )pbdoc",
         py::arg("k") = 2.)
    .def("configLambda",
         &dpd::DPDSimulator::configLambda,
         R"pbdoc(
          Configure the modified velocity verlet integration parameter `\lambda`.
     )pbdoc",
         py::arg("l") = 0.65)
    .def("configSlipspringHighCutoff",
         &dpd::DPDSimulator::configSlipspringHighCutoff,
         R"pbdoc(
          Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.
     )pbdoc",
         py::arg("cutoff") = 2.)
    .def("configSlipspringLowCutoff",
         &dpd::DPDSimulator::configSlipspringLowCutoff,
         R"pbdoc(
          Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.
     )pbdoc",
         py::arg("cutoff") = 0.5)
    .def_static("readRestartFile",
                &dpd::DPDSimulator::readRestartFile,
                R"pbdoc(
          Read a restart file in order to continue a simulation.
     )pbdoc",
                py::arg("file"))
    .def("configRestartOutput",
         &dpd::DPDSimulator::configRestartOutput,
         R"pbdoc(
          Set when to output a restart where.

          Note:
               The filename determines the type of serialisation: 
               .json, .xml are supported; other file endings will lead to binary serialisation (fastest!).

          Caution:
               This method may not be backwards- nor forward-compatible.
               Use the same version of pylimer-tools if you want to be sure that things work.

          Arguments:
               - file: the file path to the restart file to write
               - outputEvery: how often to write the restart file
     )pbdoc",
         py::arg("file"),
         py::arg("outputEvery") = 50000)
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
               - p
     )pbdoc",
         py::arg("values"),
         py::arg("numcorrin") = 32,
         py::arg("p") = 16,
         py::arg("m") = 2)
    .def("configStepOutput",
         &dpd::DPDSimulator::configStepOutput,
         R"pbdoc(
          Set which values to log.
          
          Arguments:
               - values: a list of OutputConfiguration structs
     )pbdoc")
    .def("configShiftPossibilityEmpty",
         &dpd::DPDSimulator::configShiftPossibilityEmpty,
         R"pbdoc()pbdoc")
    .def("configShiftOneAtATime",
         &dpd::DPDSimulator::configShiftOneAtATime,
         R"pbdoc()pbdoc")
    .def("configNumStepsMC", &dpd::DPDSimulator::configNumStepsMC, R"pbdoc(
          Configure the number of steps to do in one MC sequence.
     )pbdoc")
    .def("configNumStepsDPD", &dpd::DPDSimulator::configNumStepsDPD, R"pbdoc(
          Configure the number of steps to do in one DPD sequence.
     )pbdoc")
    .def("configBondFormation",
         &dpd::DPDSimulator::configBondFormation,
         R"pbdoc(
          Configure how to do bond formation during the run.

          Arguments:
          - num_bonds_to_form (int): the nr of bonds to form in total. Use 0 to stop bond formation.
          - num_bonds_per_atom_type (dict): the nr of bonds each atom type may have at most (e.g., 2 for strand atoms, 4 for a tertiary cross-links)
          - bond_formation_dist (float): the maximum distance allowed to form bonds
          - attempt_bond_formation_every (int): attempt to form bonds every this many steps during the simulation run
         )pbdoc",
         py::arg("num_bonds_to_form"),
         py::arg("max_bonds_per_atom_type"),
         py::arg("bond_formation_dist") = 1.0,
         py::arg("attempt_bond_formation_every") = 50)
     .def("getNrOfBondsToForm", &dpd::DPDSimulator::getNrOfBondsToForm, R"pbdoc(
          Get the number of bonds that are configured to have to be formed.
     )pbdoc")
    .def("configAllowRelocationInNetwork",
         &dpd::DPDSimulator::configAllowRelocationInNetwork,
         R"pbdoc(
          Configure whether a relocation step may happen when a slip-spring has ended at a cross-link.
          
          Side-effect: if true, the relocations may also happen *to* a slip-spring next to a cross-link.

          Arguments:
          - allow_relocation_in_network (bool): Whether to allow relocation in the network or not.
         )pbdoc",
         py::arg("allow_relocation_in_network") = false)
    .def("startMeasuringMSDForAtoms",
         &dpd::DPDSimulator::startMeasuringMSDForAtoms,
         R"pbdoc(
          Set a new origin for measuing the mean square displacement for a specified set of atoms
         )pbdoc",
         py::arg("atom_ids"))
    .def("getUniverse", &dpd::DPDSimulator::getUniverse, R"pbdoc(
     Get a universe instance from the current coordinates (and connectivity).

     Arguments:
          - with_slip_springs (bool): whether to include slip-springs in the returned universe.
    )pbdoc", py::arg("with_slipsprings") = true)
    .def("getTimestep", &dpd::DPDSimulator::getTimestep)
    .def("getCurrentTimestep", &dpd::DPDSimulator::getCurrentTimestep)
    .def("getTemperature", &dpd::DPDSimulator::getTemperature)
    .def("getSpringConstant", &dpd::DPDSimulator::getSpringConstant)
    .def("getShiftOneAtATime", &dpd::DPDSimulator::getShiftOneAtATime)
    .def("getNumSlipSprings", &dpd::DPDSimulator::getNumSlipSprings)
    .def("getNumStepsDPD", &dpd::DPDSimulator::getNumStepsDPD)
    .def("getNumStepsMC", &dpd::DPDSimulator::getNumStepsMC)
    .def("getShiftPossibilityEmpty",
         &dpd::DPDSimulator::getShiftPossibilityEmpty)
    .def("validateNeighbourList", &dpd::DPDSimulator::validateNeighbourlist)
    .def("validateState", &dpd::DPDSimulator::validateState);
}

#endif /* PYBIND_CALC_H */
