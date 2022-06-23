#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MEHPForceRelaxation.h"
#include "../calc/MEHPanalysis.h"
#include "../calc/MMTanalysis.h"
#include "../entities/Universe.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pylimer_tools::calc;
namespace pe = pylimer_tools::entities;

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

  py::enum_<mehp::ExitReason>(m, "ExitReason")
    .value("UNSET", mehp::ExitReason::UNSET)
    .value("MAX_STEPS", mehp::ExitReason::MAX_STEPS)
    .value("F_TOLERANCE", mehp::ExitReason::F_TOLERANCE)
    .value("X_TOLERANCE", mehp::ExitReason::X_TOLERANCE)
    .value("OTHER", mehp::ExitReason::OTHER);

  py::class_<mehp::MEHPForceRelaxation>(m, "MEHPForceRelaxation", R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
     )pbdoc")
    .def(py::init<pe::Universe, int, bool>(),
         R"pbdoc(
          Instantiate the simulator for a certain universe.

          :param universe: the universe to simulate with
          :param crosslinkerType: The atom type of the cross-linkers. Needed to reduce the network.
          :param is2D: Whether to ignore the z direction.
          )pbdoc",
         py::arg("universe"),
         py::arg("crosslinkerType") = 2,
         py::arg("is2D") = false)
    .def("runForceRelaxation",
         &mehp::MEHPForceRelaxation::runForceRelaxation,
         R"pbdoc(
          Run the simulation.

          :param algorithm: The algorithm to use for the force relaxation. Choices: see `NLopt Algorithms <https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/>`_
          :param maxNrOfSteps: The maximum number of steps to do during the simulation.
          :param xTolerance: The tolerance of the displacements as an exit condition.
          :param fTolerance: The tolerance of the force as an exit condition.
          :param loopTol: 
          :param is2d: Specify true if you want to evaluate the force relation only in x and y direction.
          )pbdoc",
         py::arg("algorithm") = "LD_MMA",
         py::arg("maxNrOfSteps") = 250000,
         py::arg("xTolerance") = 1e-12,
         py::arg("fTolerance") = 1e-9,
         py::arg("loopTol") = 1e-2,
         py::arg("kappa") = 1.0)
    .def("getForce",
         &mehp::MEHPForceRelaxation::getForce,
         R"pbdoc(
          Returns the force at the current state of the simulation.
          
          :param kappa: The spring constant to use for the force evaluation.
     )pbdoc",
         py::arg("kappa") = 1.0)
    .def("getResidualNorm",
         &mehp::MEHPForceRelaxation::getResidualNorm,
         R"pbdoc(
          Returns the residual norm at the current state of the simulation.
          
          :param kappa: The spring constant to use for the force evaluation.
     )pbdoc",
         py::arg("kappa") = 1.0)
    .def("getPressure",
         &mehp::MEHPForceRelaxation::getPressure,
         R"pbdoc(
          Returns the pressure at the current state of the simulation.
          
          :param kappa: The spring constant to use for the force evaluation.
     )pbdoc",
         py::arg("kappa") = 1.0)
    .def("getStressTensor",
         &mehp::MEHPForceRelaxation::getStressTensor,
         R"pbdoc(
          Returns the stress tensor at the current state of the simulation.
          
          :param kappa: The spring constant to use for the force evaluation.
     )pbdoc",
         py::arg("kappa") = 1.0)
    .def("getGammaFactor",
         &mehp::MEHPForceRelaxation::getGammaFactor,
         R"pbdoc(
          omputes the gamma factor as part of the ANT/MEHP formulism, i.e.:

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
     )pbdoc")
    .def("getNrOfActiveNodes",
         &mehp::MEHPForceRelaxation::getNrOfActiveNodes,
         R"pbdoc(
           Get the number of active nodes remaining after running the simulation.
     )pbdoc")
    .def("getNrOfActiveSprings",
         &mehp::MEHPForceRelaxation::getNrOfActiveSprings,
         R"pbdoc(
           Get the number of active springs remaining after running the simulation.
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
}

#endif /* PYBIND_CALC_H */
