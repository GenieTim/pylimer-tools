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

void init_pylimer_bound_calc(py::module_ &m) {
  m.def("predictGelationPoint", &mmt::predictGelationPoint,
        "Predict the gelation point of a Universe");
  // m.def("computeExtentOfReaction", &mmt::computeExtentOfReaction, "Compute
  // extent of reaction");
  m.def("computeStoichiometricInbalance", &mmt::computeStoichiometricInbalance,
        "Compute stoichiometric inbalance");

  py::class_<mehp::MEHPForceRelaxation>(m, "MEHPForceRelaxation", R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
  )pbdoc")
      .def(py::init<pe::Universe, int>())
      .def("configDoOutputSteps",
           &mehp::MEHPForceRelaxation::configDoOutputSteps, R"pbdoc(
          Whether to output progress, where, and how often.

          :param outputFile: the file to write the time-steps to.
          :param outputFrequency: how often to output the time-step. Set to 0 (default) to disable progress output.
          )pbdoc",
           py::arg("outputFile"), py::arg("outputFrequency") = 0)
      .def("configDoOutputFinalCoordinates",
           &mehp::MEHPForceRelaxation::configDoOutputFinalCoordinates)
      .def("runForceRelaxation", &mehp::MEHPForceRelaxation::runForceRelaxation,
           R"pbdoc(
          Run the simulation.

          :param crosslinkerType: The atom type of the cross-linkers. Needed to reduce the network.
          :param maxNrOfSteps: The maximum number of steps to do during the simulation.
          :param tolerance: The tolerance of the force as an exit condition.
          :param Nb2: The denominator in the equation of :math:`\Gamma` (see: :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaEq()`). If -1.0 (default), the network is used for determination.
          :param is2d: Specify true if you want to evaluate the force relation only in x and y direction.
          )pbdoc",
           py::arg("crosslinkerType") = 2, py::arg("maxNrOfSteps") = 250000,
           py::arg("tolerance") = 1e-8, py::arg("Nb2") = -1.0,
           py::arg("is2d") = false, py::arg("dt") = 0.05,
           py::arg("kappa") = 1.0)
      //  .def("getVolume", &mehp::MEHPForceRelaxation::getVolume)
      .def("getNrOfNodes", &mehp::MEHPForceRelaxation::getNrOfNodes, R"pbdoc(
           Get the number of nodes considered in this simulation.
      )pbdoc")
      .def("getNrOfSprings", &mehp::MEHPForceRelaxation::getNrOfSprings,
           R"pbdoc(
           Get the number of springs considered in this simulation.
      )pbdoc")
      .def("getNrOfActiveNodes", &mehp::MEHPForceRelaxation::getNrOfActiveNodes,
           R"pbdoc(
           Get the number of active nodes remaining after running the simulation.
      )pbdoc")
      .def("getNrOfActiveSprings",
           &mehp::MEHPForceRelaxation::getNrOfActiveSprings, R"pbdoc(
           Get the number of active springs remaining after running the simulation.
      )pbdoc")
      .def("getAverageSpringLength",
           &mehp::MEHPForceRelaxation::getAverageSpringLength, R"pbdoc(
           Get the average length of the springs.
      )pbdoc")
      .def("getGammaEq", &mehp::MEHPForceRelaxation::getGammaEq, R"pbdoc(
          Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

          :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \frac{\bar{r_{\eta}}^2}{N_{\eta} b^2}`,
          which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
          where :math:`\eta` is the index of a particular strand, 
          :math:`N_{\eta}` is the number of atoms in this strand :math:`\eta`, 
          :math:`b` its mean square bond length,
          :math:`T` the temperature and :math:`k_B` Boltzmann's constant.
      )pbdoc")
      .def("getGammaX", &mehp::MEHPForceRelaxation::getGammaX, R"pbdoc(
           Get three times the contribution of the x-direction to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaEq()`.
      )pbdoc")
      .def("getGammaY", &mehp::MEHPForceRelaxation::getGammaY, R"pbdoc(
           Get three times the contribution of the y-direction to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaEq()`.
      )pbdoc")
      .def("getGammaZ", &mehp::MEHPForceRelaxation::getGammaZ, R"pbdoc(
           Get three times the contribution of the z-direction to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaEq()`.
      )pbdoc")
      .def("getNb2", &mehp::MEHPForceRelaxation::getNb2, R"pbdoc(
           Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaEq()` for :math:`N_{\eta} b^2`.
      )pbdoc");
}

#endif /* PYBIND_CALC_H */
