#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MEHPForceBalance.h"
#include "../calc/MEHPForceBalance2.h"
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
    .def(py::init<pe::Universe, int, bool, mehp::MEHPForceEvaluator*>(),
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
         py::arg("forceEvaluator") = nullptr)
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
    .def("getEffectiveFunctionalityOfAtoms",
         &mehp::MEHPForceRelaxation::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
          Returns the number of active springs connected to each atom, atomId used as index

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
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
     )pbdoc",
         py::arg("newCrosslinkerType") = 2);

  py::enum_<mehp::BalanceRunMode>(m, "BalanceRunMode")
    .value("EIGEN_ALL", mehp::BalanceRunMode::EIGEN_ALL)
    .value("EIGEN_RANDOM", mehp::BalanceRunMode::EIGEN_RANDOM)
    .value("EIGEN_HEURISTIC", mehp::BalanceRunMode::EIGEN_HEURISTIC)
    .value("EIGEN_STRANDS", mehp::BalanceRunMode::EIGEN_STRANDS)
    .value("ITERATIVE", mehp::BalanceRunMode::ITERATIVE);

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
         py::arg("remove2functionalCrosslinkers") = false)
    .def_property_readonly("network", &mehp::MEHPForceBalance::getNetwork)
    .def("validateNetwork",
         py::overload_cast<>(&mehp::MEHPForceBalance::validateNetwork),
         R"pbdoc(
           Validates the internal structures.
     )pbdoc")
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
         py::arg("innerMaxNrOfSteps") = 250,
         py::arg("innerAlphaTolerance") = 1e-15,
         py::arg("oneOverSpringPartitionUpperLimit") = 1.0,
         py::arg("maxFlag") = 7,
         py::arg("allowRemovalOfSliplinks") = false,
         py::arg("allowMoveOfSliplinks") = false)
    .def("inspectLinkDisplacementToMeanPositionUpdate",
         &mehp::MEHPForceBalance::inspectLinkDisplacementToMeanPositionUpdate,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("damping") = 1.0)
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
         py::arg("excludeCrosslinks") = false)
    .def("getStressTensor",
         &mehp::MEHPForceBalance::getStressTensor,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc",
         py::arg("oneOverSpringPartitionUpperLimit") = 1.)
    .def("getStressTensorLinkBased",
         &mehp::MEHPForceBalance::getStressTensorLinkBased,
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
     )pbdoc",
         py::arg("newCrosslinkerType") = 2);

  ////////////////////////////////////////////////////////////////

  py::class_<mehp::MEHPForceBalance2>(m, "MEHPForceBalance2", R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
     )pbdoc")
    .def(py::init<pe::Universe, int, bool>(),
         R"pbdoc(
          Instantiate the simulator for a certain universe.

          :param universe: the universe to simulate with
          :param crosslinkerType: The atom type of the cross-linkers. Needed to reduce the network.
          :param is2D: Whether to ignore the z direction.
          :param forceEvaluator: The force evaluator to use
          )pbdoc",
         py::arg("universe"),
         py::arg("crosslinkerType") = 2,
         py::arg("is2D") = false)
    .def_property_readonly("network", &mehp::MEHPForceBalance2::getNetwork)
    .def("validateNetwork",
         py::overload_cast<>(&mehp::MEHPForceBalance2::validateNetwork),
         R"pbdoc(
           Validates the internal structures.
     )pbdoc")
    .def("runForceRelaxation",
         &mehp::MEHPForceBalance2::runForceRelaxation,
         R"pbdoc(
          Run the simulation.
          Note that the final state of the minimization is persisted and reused if you use this method again.
          This is useful if you want to run a global optimization first and add a local one afterwards.
          As a consequence though, you cannot simply benchmark only this method; you must include the setup.
          )pbdoc",
         py::arg("algorithm") = "LD_MMA",
         py::arg("maxNrOfSteps") = 250000,
         py::arg("xTolerance") = 1e-12,
         py::arg("fTolerance") = 1e-9,
         py::arg("constraintTol") = 1e-9)
    .def("inspectLinkDisplacementToMeanPositionUpdate",
         &mehp::MEHPForceBalance2::inspectLinkDisplacementToMeanPositionUpdate,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("damping") = 1.0)
    .def("getForceOn",
         &mehp::MEHPForceBalance2::getForceOn,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
    .def("inspectDisplacementToMeanPositionUpdate",
         &mehp::MEHPForceBalance2::inspectDisplacementToMeanPositionUpdate,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
    .def("inspectSpringPartitionUpdate",
         &mehp::MEHPForceBalance2::inspectSpringPartitionUpdate,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"))
    .def("inspectParametrisationOptimsationForLink",
         &mehp::MEHPForceBalance2::inspectParametrisationOptimsationForLink,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"),
         py::arg("displacements"),
         py::arg("springPartitions"),
         py::arg("maxNrOfSteps") = 100,
         py::arg("alpha_tol") = 1e-9,
         py::arg("distanceBackTolerance") = 1e-9,
         py::arg("residualNormSTolerance") = 1e-20,
         py::arg("minNrOfSteps") = 1,
         py::arg("oneOverSpringPartitionUpperLimit") = 1.0)
    .def("getSpringpartitionIndicesOfSliplink",
         &mehp::MEHPForceBalance2::getSpringpartitionIndicesOfSliplink,
         R"pbdoc()pbdoc",
         py::arg("linkIdx"))
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
    .def("getPressure",
         &mehp::MEHPForceBalance2::getPressure,
         R"pbdoc(
          Returns the pressure at the current state of the simulation.
     )pbdoc")
    .def("addSlipLinks",
         py::overload_cast<const std::vector<size_t>,
                           const std::vector<size_t>,
                           const std::vector<double>,
                           const std::vector<double>,
                           const std::vector<double>,
                           const std::vector<double>,
                           const std::vector<double>,
                           const bool>(&mehp::MEHPForceBalance2::addSlipLinks),
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
    .def("getStressTensor",
         &mehp::MEHPForceBalance2::getStressTensor,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc",
         py::arg("oneOverSpringPartitionUpperLimit") = -1)
    .def("getStressTensorLinkBased",
         &mehp::MEHPForceBalance2::getStressTensorLinkBased,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
     )pbdoc",
         py::arg("oneOverSpringPartitionUpperLimit") = -1,
         py::arg("xlinksOnly") = false)
    .def("getGammaFactor",
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
          
          :param r0squared: The denominator in the equation of :math:`\Gamma`. If :math:`-1.0` (default), the network is used for determination (which is not accurate). For phantom systems, the correct value is :math:`Nb^2`.
               For other systems, the value could be determined by `~pylimer_tools_cpp.pylimer_tools_cpp.Universe.computeMeanEndToEndDistance` on the melt system.
          :param nrOfChains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains. 
     )pbdoc",
         py::arg("r0squared") = -1.0,
         py::arg("nrOfChains") = -1)
    .def("getNrOfNodes", &mehp::MEHPForceBalance2::getNrOfNodes, R"pbdoc(
           Get the number of nodes considered in this simulation.
     )pbdoc")
    .def("getNrOfSprings",
         &mehp::MEHPForceBalance2::getNrOfSprings,
         R"pbdoc(
          Get the number of springs considered in this simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc")
    .def("getSpringPartitions",
         &mehp::MEHPForceBalance2::getSpringPartitions,
         R"pbdoc(
          Get the current spring partitions.
     )pbdoc")
    .def("setSpringPartitions",
         &mehp::MEHPForceBalance2::setSpringPartitions,
         R"pbdoc(
          Set the current spring partitions.
     )pbdoc")
    .def("getDisplacements",
         &mehp::MEHPForceBalance2::getCurrentDisplacements,
         R"pbdoc(
          Get the current link displacements.
     )pbdoc")
    .def("setDisplacements",
         &mehp::MEHPForceBalance2::setCurrentDisplacements,
         R"pbdoc(
          Set the current link displacements.
     )pbdoc")
    .def("setSpringContourLengths",
         &mehp::MEHPForceBalance2::setSpringContourLengths,
         R"pbdoc(
          Set/overwrite the contour lengths.
     )pbdoc")
    .def("getDisplacementResidualNorm",
         py::overload_cast<double>(
           &mehp::MEHPForceBalance2::getDisplacementResidualNorm),
         R"pbdoc(
          Get the current link displacement residual norm.
     )pbdoc")
    .def("getIdsOfActiveNodes",
         &mehp::MEHPForceBalance2::getIdsOfActiveNodes,
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
         &mehp::MEHPForceBalance2::getNrOfActiveNodes,
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
         &mehp::MEHPForceBalance2::getNrOfActiveSprings,
         R"pbdoc(
           Get the number of active springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getNrOfActivePartialSprings",
         &mehp::MEHPForceBalance2::getNrOfActivePartialSprings,
         R"pbdoc(
           Get the number of active partial springs remaining after running the simulation.

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getEffectiveFunctionalityOfAtoms",
         &mehp::MEHPForceBalance2::getEffectiveFunctionalityOfAtoms,
         R"pbdoc(
          Returns the number of active springs connected to each atom, atomId used as index

          :param tolerance: springs under this length are considered inactive
     )pbdoc",
         py::arg("tolerance") = 0.1)
    .def("getAverageSpringLength",
         &mehp::MEHPForceBalance2::getAverageSpringLength,
         R"pbdoc(
           Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()`,
           this value is normalized by the number of springs rather than the number of chains.
     )pbdoc")
    .def("getDefaultR0Square",
         &mehp::MEHPForceBalance2::getDefaultR0Square,
         R"pbdoc(
           Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()` for :math:`\langle R_{0,\eta}^2\rangle`.
     )pbdoc")
    .def("getDefaultNrOfChains",
         &mehp::MEHPForceBalance2::getDefaultNrOfChains,
         R"pbdoc(
          Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()` for normalizing the distances.`.
     )pbdoc")
    .def("getNrOfIterations",
         &mehp::MEHPForceBalance2::getNrOfIterations,
         R"pbdoc(
          Returns the number of iterations used for force relaxation.
     )pbdoc")
    .def("getExitReason", &mehp::MEHPForceBalance2::getExitReason, R"pbdoc(
           Returns the reason for termination of the simulation
     )pbdoc")
    .def("getCrosslinkerVerse",
         &mehp::MEHPForceBalance2::getCrosslinkerVerse,
         R"pbdoc(
          Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
     )pbdoc",
         py::arg("newCrosslinkerType") = 2);
}

#endif /* PYBIND_CALC_H */
