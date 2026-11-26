import copy
import itertools
import os
import sys
import unittest

import numpy as np
from pint import UnitRegistry

from pylimer_tools.calc.miller_macosko_theory import (
    _compute_validate_parameters,
    compute_entanglement_modulus,
    compute_miller_macosko_probabilities,
    compute_modulus_decomposition,
    compute_weight_fraction_of_soluble_material,
    compute_weight_fraction_of_soluble_material_from_weight_fractions,
    predict_gelation_point,
    predict_number_density_of_junction_points,
    predict_p_from_w_sol,
    compute_weight_fraction_of_backbone,
    compute_weight_fraction_of_dangling_chains,
    predict_number_density_of_network_strands,
    predict_shear_modulus,
    compute_extracted_modulus,
    compute_junction_modulus,
    predict_maximum_p,
)
from pylimer_tools.calc.structure_analysis import (
    compute_extent_of_reaction,
    compute_stoichiometric_imbalance,
    compute_weight_fractions,
    measure_weight_fraction_of_backbone,
    measure_weight_fraction_of_dangling_chains,
    measure_lower_bound_weight_fraction_of_soluble_material,
    measure_weight_fraction_of_soluble_material,
)
from pylimer_tools.io.unit_styles import UnitStyleFactory

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))
from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMMTAnalysisFunctions(UniverseUsingTestCase):
    def test_stoichiometric_imbalance(self):
        self.assertAlmostEqual(
            0, compute_stoichiometric_imbalance(self.emptyUniverse, 2)
        )
        self.assertAlmostEqual(
            (3 * 3) / (3 * 2),
            compute_stoichiometric_imbalance(
                self.testUniverse, 2, functionality_per_type={1: 2, 2: 3}
            ),
        )
        self.assertAlmostEqual(
            (3 * 3) / (3 * 2),
            compute_stoichiometric_imbalance(
                self.testUniverse,
                2,
            ),
        )
        self.assertAlmostEqual(
            (4 * 3) / (3 * 2),
            compute_stoichiometric_imbalance(
                self.testUniverse, 2, functionality_per_type={1: 2, 2: 4}
            ),
        )
        self.assertAlmostEqual(
            (1 * 2 + 1 * 3 + 0) / (3 * 2),
            compute_stoichiometric_imbalance(self.testUniverse, 2, effective=True),
        )
        self.assertAlmostEqual(
            0.0, compute_stoichiometric_imbalance(self.testUniverse, 7)
        )
        self.assertAlmostEqual(
            compute_stoichiometric_imbalance(
                self.testUniverse, 2, functionality_per_type={1: 2, 2: 3}
            ),
            ((3 * 3) / ((5 * 2) / (5 / 3))),
        )
        self.testUniverse.add_atoms([100], [3], [0], [0], [0], [0], [0], [0])
        self.assertAlmostEqual(
            compute_stoichiometric_imbalance(self.testUniverse, 2, ignore_types=[3]),
            ((3 * 3) / ((5 * 2) / (5 / 3))),
        )

    def test_extent_of_reaction(self):
        self.assertAlmostEqual(1.0, compute_extent_of_reaction(self.emptyUniverse, 2))
        self.assertAlmostEqual(
            5.0 / 6.0, compute_extent_of_reaction(self.testUniverse, 2)
        )

    def test_gelation_point_prediction(self):
        self.assertAlmostEqual(1, predict_gelation_point(1, 2))
        self.assertAlmostEqual(1, predict_gelation_point(1, 2, 1))

    def test_shear_modulus_prediction(self):
        self.assertRaises(
            ValueError,
            lambda: predict_shear_modulus(
                network=self.emptyUniverse, crosslinker_type=2, unit_style=None
            ),
        )
        self.assertRaises(
            ValueError,
            lambda: predict_shear_modulus(network=self.emptyUniverse, unit_style=None),
        )
        self.saturatedTestUniverse.set_masses({1: 1, 2: 1})

        unit_style_factory = UnitStyleFactory()
        unit_style = unit_style_factory.get_unit_style("si")
        # TODO: find literature motiviation for results fo the functions
        self.assertAlmostEqual(
            0.1467020757993193,
            predict_shear_modulus(
                network=self.saturatedTestUniverse,
                unit_style=unit_style,
                crosslinker_type=2,
            )
            .to("MPa")
            .magnitude,
        )
        self.assertAlmostEqual(
            0.12832791725960474,
            predict_shear_modulus(
                network=self.saturatedTestUniverse,
                unit_style=unit_style,
                crosslinker_type=2,
                functionality_per_type={2: 4},
            )
            .to("MPa")
            .magnitude,
        )
        self.saturatedTestUniverse.set_masses({1: 1, 2: 1})
        self.assertAlmostEqual(
            0.1467020757993193,
            predict_shear_modulus(
                network=self.saturatedTestUniverse,
                unit_style=unit_style,
                crosslinker_type=2,
            )
            .to("MPa")
            .magnitude,
        )
        ureg = UnitRegistry()
        self.assertAlmostEqual(
            compute_junction_modulus(
                p=1, r=1, f=4, xlink_concentration_0=0 * ureg("g/nm^3"), b2=1, ureg=ureg
            ).magnitude,
            0.0,
        )
        self.assertAlmostEqual(
            compute_extracted_modulus(
                p=1,
                r=1,
                f=4,
                b2=1,
                xlink_concentration_0=0 * ureg("1/nm^3"),
                w_sol=0.1,
                ureg=ureg,
                g_e_1=0.3 * ureg("MPa"),
            ).magnitude,
            0.0,
        )

    def test_predict_number_density_of_junction_points(self):
        self.testUniverse.set_masses({1: 1, 2: 1})
        self.assertAlmostEqual(
            1.5, compute_stoichiometric_imbalance(self.testUniverse, 2)
        )
        # TODO: find literature motiviation for results fo the functions
        self.assertAlmostEqual(
            predict_number_density_of_junction_points(self.testUniverse, 2), 0.0
        )
        self.assertRaises(
            NotImplementedError,
            lambda: predict_number_density_of_junction_points(
                self.testUniverse, 2, functionality_per_type={1: 2, 2: 5}
            ),
        )
        self.assertAlmostEqual(
            predict_number_density_of_junction_points(
                self.testUniverse, 2, functionality_per_type={1: 2, 2: 3}
            ),
            0.0,
        )
        self.saturatedTestUniverse.set_masses({1: 1, 2: 1})
        self.assertAlmostEqual(
            predict_number_density_of_junction_points(
                self.saturatedTestUniverse, 2, functionality_per_type={1: 2, 2: 4}
            ),
            0.08304616298035744,
        )

    def test_predict_number_density_of_network_strands(self):
        self.testUniverse.set_masses({1: 1, 2: 1})
        # TODO: find literature motiviation for results fo the functions
        self.assertAlmostEqual(
            predict_number_density_of_network_strands(self.testUniverse, 2), 0.0
        )
        self.assertRaises(
            NotImplementedError,
            lambda: predict_number_density_of_network_strands(
                self.testUniverse, 2, functionality_per_type={1: 2, 2: 5}
            ),
        )
        self.assertAlmostEqual(
            predict_number_density_of_network_strands(
                self.testUniverse, 2, functionality_per_type={1: 2, 2: 3}
            ),
            0.0,
        )
        self.saturatedTestUniverse.set_masses({1: 1, 2: 1})
        self.assertAlmostEqual(
            predict_number_density_of_network_strands(
                self.saturatedTestUniverse, 2, functionality_per_type={1: 2, 2: 4}
            ),
            0.13752853164771697,
        )

    def test_weight_fraction_calculations(self):
        self.assertDictEqual({}, compute_weight_fractions(self.emptyUniverse))

        self.testUniverse.set_masses({1: 1, 2: 1})
        weight_fractions = compute_weight_fractions(self.testUniverse)
        self.assertDictEqual(weight_fractions, {1: 1 - 3.0 / 8.0, 2: 3.0 / 8.0})
        test_universe_copy = copy.copy(self.testUniverse)
        test_universe_copy.remove_atoms([1, 2, 3, 4, 5, 6])
        self.assertTrue(test_universe_copy.get_nr_of_atoms() == 2)
        self.assertDictEqual(
            compute_weight_fractions(test_universe_copy), {1: 1.0 / 2.0, 2: 1.0 / 2.0}
        )

        params = {
            "crosslinker_type": 2,
            "functionality_per_type": {1: 2, 2: 3},
            "r": 1.0,
            "p": 0.95,
            "weight_fractions": {1: 0.85, 2: 0.15},
            "b2": 1,
        }

        ps = [1.0, 0.9, 0.95, 0.98]
        rs = [1.0, 0.9, 1.1]
        fs = [3, 4, 5]

        for p, r, f in itertools.product(ps, rs, fs):
            params["p"] = min(p, predict_maximum_p(r, f))
            params["r"] = r
            params["functionality_per_type"][2] = f
            self.assertAlmostEqual(
                (
                    compute_weight_fraction_of_backbone(None, **params)
                    + compute_weight_fraction_of_dangling_chains(None, **params)
                    + compute_weight_fraction_of_soluble_material(None, **params)
                ),
                1.0,
            )

    def test_max_p_prediction(self):
        for r in np.linspace(0.5, 1.0, 10):
            for f in range(3, 12):
                self.assertAlmostEqual(predict_maximum_p(r=r, f=f, b2=1), 1.0)
        for r in np.linspace(1.01, 1.5, 10):
            for f in range(3, 12):
                self.assertLess(predict_maximum_p(r=r, f=f, b2=1), 1.0)

    def test_soluble_weight_fraction_measurement(self):
        self.testUniverse.set_masses({1: 1, 2: 1})
        self.assertEqual(
            measure_weight_fraction_of_soluble_material(self.emptyUniverse), 0.0
        )
        self.assertEqual(
            measure_weight_fraction_of_soluble_material(self.testUniverse, rel_tol=0),
            0.0,
        )
        self.assertEqual(
            measure_weight_fraction_of_soluble_material(
                self.testUniverse, abs_tol=10000
            ),
            1.0,
        )
        self.assertEqual(
            measure_weight_fraction_of_soluble_material(self.testUniverse), 1 / 8
        )
        self.assertEqual(
            measure_lower_bound_weight_fraction_of_soluble_material(
                self.emptyUniverse, 2
            ),
            0.0,
        )
        self.assertEqual(
            measure_lower_bound_weight_fraction_of_soluble_material(
                self.testUniverse, 2, abs_tol=0
            ),
            0.0,
        )
        self.assertEqual(
            measure_lower_bound_weight_fraction_of_soluble_material(
                self.testUniverse, 2, abs_tol=10000
            ),
            0.125,
        )

    def test_soluble_material_weight_fraction_calculation(self):
        self.testUniverse.set_masses({1: 1, 2: 1})
        self.assertRaises(
            NotImplementedError,
            lambda: compute_weight_fraction_of_soluble_material(
                self.testUniverse, 2, functionality_per_type={1: 2, 2: 2}
            ),
        )
        self.assertRaises(
            ValueError,
            lambda: compute_weight_fraction_of_soluble_material(self.testUniverse, 7),
        )
        self.saturatedTestUniverse.set_masses({1: 1, 2: 1})

        self.assertAlmostEqual(
            0.010699588477366243,
            compute_weight_fraction_of_soluble_material(self.saturatedTestUniverse, 2),
        )
        self.assertAlmostEqual(
            predict_p_from_w_sol(
                w_sol=0,
                network=None,
                crosslinker_type=2,
                functionality_per_type={1: 2, 2: 4},
                weight_fractions={1: 1.0, 2: 0.0},
                r=1,
                b2=1.0,
            ),
            1.0,
            3,
        )
        # Verify compute_weight_fraction_of_soluble_material_from_weight_fractions
        w_sol_dict = compute_weight_fraction_of_soluble_material(
            network=None,
            crosslinker_type=2,
            functionality_per_type={1: 2, 2: 4},
            weight_fractions={1: 0.85, 2: 0.15},
            r=1.0,
            p=0.95,
            b2=1.0,
        )
        w_sol_direct = (
            compute_weight_fraction_of_soluble_material_from_weight_fractions(
                r=1.0, p=0.95, f=4, w_f=0.15, w_g=0.85, g=2
            )
        )
        self.assertAlmostEqual(w_sol_dict, w_sol_direct, places=10)

        # Verify parameter validation error handling in _compute_validate_parameters
        with self.assertRaises(ValueError):
            _compute_validate_parameters({"r": -0.5, "p": 0.9}, ["r", "p"])
        with self.assertRaises(ValueError):
            _compute_validate_parameters({"r": 1.0, "p": 1.5}, ["r", "p"])
        with self.assertRaises(ValueError):
            _compute_validate_parameters({}, ["r", "p"])

    # def testProbabilityCalculations(self):
    #     self.assertRaises(
    #         ValueError, lambda: compute_miller_macosko_probabilities(0.9, 2, 2))
    #     self.assertRaises(
    # ValueError, lambda: compute_miller_macosko_probabilities(0.1, 0.9, 3))

    def test_backbone_weight_fraction_calculations(self):
        self.assertEqual(0, compute_weight_fraction_of_backbone(self.emptyUniverse, 2))
        self.assertEqual(
            0.0, compute_weight_fraction_of_dangling_chains(self.emptyUniverse, 2)
        )
        self.assertEqual(
            1, compute_weight_fraction_of_soluble_material(self.emptyUniverse, 2)
        )
        self.saturatedTestUniverse.set_masses({1: 1, 2: 0})
        self.assertAlmostEqual(
            1.0, compute_extent_of_reaction(self.saturatedTestUniverse, 2)
        )
        # TODO: get some literature backed values to test for
        bb = compute_weight_fraction_of_backbone(
            self.saturatedTestUniverse, crosslinker_type=2
        )
        wsol = compute_weight_fraction_of_soluble_material(
            self.saturatedTestUniverse, crosslinker_type=2
        )
        self.assertAlmostEqual(0.8, bb, places=1)
        self.saturatedTestUniverse.set_masses({1: 1, 2: 0})
        self.assertAlmostEqual(
            1 - bb - wsol,
            compute_weight_fraction_of_dangling_chains(self.saturatedTestUniverse, 2),
            places=5,
        )

        self.saturatedTestUniverse.set_masses({1: 1, 2: 1})
        # test also as if the functionality was 4
        # self.assertRaises(ValueError, lambda: compute_weight_fraction_of_backbone(self.saturatedTestUniverse,
        #       crosslinker_type=2, functionality_per_type={
        #     1: 2, 2: 4
        # }))
        # NOTE: requires a short strand length with these systems, as otherwise, r > 1
        # which is not supported by the formulas implemented
        self.assertAlmostEqual(
            0.7205133203806013,
            compute_weight_fraction_of_backbone(
                self.saturatedTestUniverse,
                crosslinker_type=2,
                functionality_per_type={1: 2, 2: 4},
            ),
        )

    def test_modulus_decompositions(self):
        unit_style_factory = UnitStyleFactory()
        unit_style = unit_style_factory.get_unit_style("si")
        self.assertAlmostEqual(unit_style.kb.to("J/K").magnitude, 1.381e-23)

        # these results are pretty certain, align with experimental results,
        # confirmed
        g_mmt_phantom, g_mmt_entanglement, g_anm, g_pnm = compute_modulus_decomposition(
            network=None,
            unit_style=unit_style,
            crosslinker_type=2,
            r=1.0,
            p=0.95,
            f=4,
            nu=4.69218e25 * (unit_style.get_underlying_unit_registry()("meter") ** -3),
            temperature=298 * unit_style.get_underlying_unit_registry()("kelvin"),
        )

        alpha, beta = compute_miller_macosko_probabilities(r=1.0, p=0.95, f=4)
        self.assertAlmostEqual(alpha, 0.0983588, places=5)

        self.assertAlmostEqual(g_anm.to("MPa").magnitude, 0.1930520, places=5)
        self.assertAlmostEqual(g_pnm.to("MPa").magnitude, 0.0965260, places=5)
        self.assertAlmostEqual(
            g_mmt_entanglement.to("MPa").magnitude, 0.190576, places=5
        )
        self.assertAlmostEqual(g_mmt_phantom.to("MPa").magnitude, 0.0777123, places=5)

        # these in turn require further investigation into the involvement of r
        g_mmt_phantom, g_mmt_entanglement, g_anm, g_pnm = compute_modulus_decomposition(
            network=None,
            unit_style=unit_style,
            crosslinker_type=2,
            r=1.3,
            p=0.6465,
            f=4,
            nu=1.25981e25 * (unit_style.get_underlying_unit_registry()("meter") ** -3),
            temperature=298 * unit_style.get_underlying_unit_registry()("kelvin"),
        )

        self.assertAlmostEqual(g_anm.to("MPa").magnitude, 0.051833, places=5)
        self.assertAlmostEqual(g_pnm.to("MPa").magnitude, 0.025916, places=5)
        self.assertAlmostEqual(
            g_mmt_entanglement.to("MPa").magnitude, 0.05801087, places=5
        )
        self.assertAlmostEqual(g_mmt_phantom.to("MPa").magnitude, 0.00492674, places=5)

    def test_probability_calculations(self):
        # most already in other tests
        alpha, beta = compute_miller_macosko_probabilities(r=1.0, p=0.85, f=5)
        self.assertAlmostEqual(0.282074, alpha, places=5)
        alpha, beta = compute_miller_macosko_probabilities(r=1.0, p=0.85, f=6)
        self.assertAlmostEqual(0.278715, alpha, places=5)
        self.assertRaises(
            ValueError, lambda: compute_miller_macosko_probabilities(r=1.0, p=2, f=7)
        )

    def test_monofunctional_weight_fractions_and_mass_conservation(self):
        """Test MMT weight fractions with monofunctional chains (b2 < 1) and solvent."""
        for b2 in [0.5, 0.8]:
            params = {
                "crosslinker_type": 2,
                "functionality_per_type": {1: 2, 2: 4, 3: 1},
                "weight_fractions": {1: 0.70, 2: 0.15, 3: 0.15},
                "r": 1.0,
                "p": 0.95,
                "b2": b2,
            }
            w_bb = compute_weight_fraction_of_backbone(None, **params)
            w_dang = compute_weight_fraction_of_dangling_chains(None, **params)
            w_sol = compute_weight_fraction_of_soluble_material(None, **params)
            self.assertAlmostEqual(w_bb + w_dang + w_sol, 1.0, places=10)

            # Compare backbone with same system where monofunctional has weight 0
            # (f=1 must contribute 0.0 to backbone)
            params_no_mono = {
                **params,
                "weight_fractions": {1: 0.70, 2: 0.15, 3: 0.0},
            }
            self.assertAlmostEqual(
                w_bb,
                compute_weight_fraction_of_backbone(None, **params_no_mono),
                places=10,
            )

        # Solvent (f=0) must also contribute 0.0 to backbone
        params_with_solvent = {
            "crosslinker_type": 2,
            "functionality_per_type": {1: 2, 2: 4, 4: 0},
            "weight_fractions": {1: 0.70, 2: 0.15, 4: 0.15},
            "r": 1.0,
            "p": 0.95,
            "b2": 1.0,
        }
        params_no_solvent = {
            **params_with_solvent,
            "weight_fractions": {1: 0.70, 2: 0.15, 4: 0.0},
        }
        self.assertAlmostEqual(
            compute_weight_fraction_of_backbone(None, **params_with_solvent),
            compute_weight_fraction_of_backbone(None, **params_no_solvent),
            places=10,
        )

    def test_compute_entanglement_modulus_with_beta(self):
        """Test compute_entanglement_modulus with explicit beta and without temperature."""
        ureg = UnitRegistry()
        g_e_1 = 0.5 * ureg("MPa")
        beta = 0.2
        result = compute_entanglement_modulus(g_e_1, beta=beta)
        self.assertAlmostEqual(
            result.to("MPa").magnitude, 0.5 * (1 - beta) ** 4, places=10
        )

        # Also works when computing beta from p, r, f, b2
        result_from_params = compute_entanglement_modulus(
            g_e_1, p=0.95, r=1.0, f=4, b2=1.0
        )
        _, expected_beta = compute_miller_macosko_probabilities(r=1.0, p=0.95, f=4)
        self.assertAlmostEqual(
            result_from_params.to("MPa").magnitude,
            0.5 * (1 - expected_beta) ** 4,
            places=10,
        )

        with self.assertRaises(ValueError):
            compute_entanglement_modulus(g_e_1, beta=1.5)

    def test_predict_p_from_w_sol_edge_cases(self):
        """Test predict_p_from_w_sol with r > 1 and with a 0-bond network."""
        from pylimer_tools_cpp import Universe

        # r > 1: p_max < 1.0
        r = 1.3
        f = 4
        p_true = 0.7
        weight_fractions = {1: 0.85, 2: 0.15}
        func_per_type = {1: 2, 2: f}
        w_sol = compute_weight_fraction_of_soluble_material(
            network=None,
            crosslinker_type=2,
            functionality_per_type=func_per_type,
            weight_fractions=weight_fractions,
            r=r,
            p=p_true,
            b2=1.0,
        )
        p_pred = predict_p_from_w_sol(
            w_sol=w_sol,
            network=None,
            crosslinker_type=2,
            functionality_per_type=func_per_type,
            weight_fractions=weight_fractions,
            r=r,
            b2=1.0,
        )
        self.assertAlmostEqual(p_pred, p_true, places=5)

        # 0-bond network should not short-circuit compute_wsol to 1.0
        zero_bond_net = Universe(10, 10, 10)
        zero_bond_net.add_atoms(
            [1, 2], [1, 2], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0, 0], [0, 0], [0, 0]
        )
        p_pred_net = predict_p_from_w_sol(
            w_sol=w_sol,
            network=zero_bond_net,
            crosslinker_type=2,
            functionality_per_type=func_per_type,
            weight_fractions=weight_fractions,
            r=r,
            b2=1.0,
        )
        self.assertAlmostEqual(p_pred_net, p_true, places=5)

    def test_predict_number_density_with_explicit_crosslinker_density_and_b2(self):
        """Test crosslinker_number_density and b2 forwarding in number density predictors."""
        mu = predict_number_density_of_junction_points(
            network=self.saturatedTestUniverse,
            crosslinker_type=2,
            functionality_per_type={1: 2, 2: 4},
            b2=0.9,
            crosslinker_number_density=2.5,
        )
        nu = predict_number_density_of_network_strands(
            network=self.saturatedTestUniverse,
            crosslinker_type=2,
            functionality_per_type={1: 2, 2: 4},
            b2=0.9,
            crosslinker_number_density=2.5,
        )
        self.assertGreater(mu, 0.0)
        self.assertGreater(nu, mu)

    def test_predict_number_density_weight_fractions_fallback_deprecated(self):
        """Omitting crosslinker_number_density emits DeprecationWarning."""
        self.saturatedTestUniverse.set_masses({1: 1, 2: 1})
        with self.assertWarns(DeprecationWarning):
            predict_number_density_of_junction_points(
                self.saturatedTestUniverse, 2, functionality_per_type={1: 2, 2: 4}
            )
        with self.assertWarns(DeprecationWarning):
            predict_number_density_of_network_strands(
                self.saturatedTestUniverse, 2, functionality_per_type={1: 2, 2: 4}
            )

    def test_modulus_decomposition_matches_junction_modulus_with_b2(self):
        """compute_modulus_decomposition must agree with compute_junction_modulus for b2 < 1."""
        from pylimer_tools_cpp import Universe

        unit_style = UnitStyleFactory().get_unit_style("si")
        ureg = unit_style.get_underlying_unit_registry()
        temperature = 298.15 * ureg("kelvin")
        nu = 6.0e25 * (ureg("meter") ** -3)

        for f in [3, 4]:
            for r in [0.9, 1.0, 1.1]:
                for b2 in [0.6, 0.8, 0.9, 1.0]:
                    p = min(0.95, predict_maximum_p(r=r, f=f, b2=b2))
                    # [A_f]_0 = (2 * r / (f * (2 - b2))) * nu
                    xlink_conc_0 = (2 * r / (f * (2 - b2))) * nu
                    g_mmt_phantom, _, _, _ = compute_modulus_decomposition(
                        network=None,
                        unit_style=unit_style,
                        crosslinker_type=2,
                        r=r,
                        p=p,
                        f=f,
                        nu=nu,
                        temperature=temperature,
                        b2=b2,
                    )
                    g_junction = compute_junction_modulus(
                        p=p,
                        r=r,
                        xlink_concentration_0=xlink_conc_0,
                        ureg=ureg,
                        f=f,
                        temperature=temperature,
                        b2=b2,
                    )
                    self.assertAlmostEqual(
                        g_mmt_phantom.to("MPa").magnitude,
                        g_junction.to("MPa").magnitude,
                        places=10,
                        msg=f"Mismatch for f={f}, r={r}, b2={b2}",
                    )

        # Also verify end-to-end on a Universe with B2, multi-bead B1, and multi-bead B0
        net = Universe(10, 10, 10)
        # 2 crosslinkers (type 2, f=4): ids 1, 2
        # 3 B2 chains (type 1, 3 beads each): ids 10..18 -> 6 B sites
        # 2 B1 chains (type 3, 3 beads each): ids 20..25 -> 2 B sites
        # 2 B0 solvent chains (type 4, 3 beads each): ids 30..35 -> 0 B sites
        # Total B sites = 8, total A sites = 8 => r = 1.0, b2 = 6/8 = 0.75, p = 1.0
        atom_ids = (
            [1, 2] + list(range(10, 19)) + list(range(20, 26)) + list(range(30, 36))
        )
        atom_types = [2, 2] + [1] * 9 + [3] * 6 + [4] * 6
        n_atoms = len(atom_ids)
        net.add_atoms(
            atom_ids,
            atom_types,
            [float(i) for i in range(n_atoms)],
            [0.0] * n_atoms,
            [0.0] * n_atoms,
            [0] * n_atoms,
            [0] * n_atoms,
            [0] * n_atoms,
        )
        # Internal chain bonds
        bonds_from = [10, 11, 13, 14, 16, 17, 20, 21, 23, 24, 30, 31, 33, 34]
        bonds_to = [11, 12, 14, 15, 17, 18, 21, 22, 24, 25, 31, 32, 34, 35]
        # Crosslinker bonds:
        # B2 chains (10-11-12, 13-14-15, 16-17-18) all connect xlink 1 and xlink 2
        bonds_from += [1, 12, 1, 15, 1, 18]
        bonds_to += [10, 2, 13, 2, 16, 2]
        # B1 chains (20-21-22 -> xlink 1; 23-24-25 -> xlink 2)
        bonds_from += [1, 2]
        bonds_to += [20, 23]
        net.add_bonds(bonds_from, bonds_to)
        net.set_masses({1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0})

        g_mmt_phantom, _, _, _ = compute_modulus_decomposition(
            network=net,
            unit_style=unit_style,
            crosslinker_type=2,
            temperature=temperature,
        )
        xlink_conc_net = 2 / (net.get_volume() * unit_style.get_base_unit_of("volume"))
        g_junction_net = compute_junction_modulus(
            p=1.0,
            r=1.0,
            xlink_concentration_0=xlink_conc_net,
            ureg=ureg,
            f=4,
            temperature=temperature,
            b2=0.75,
        )
        self.assertAlmostEqual(
            g_mmt_phantom.to("MPa").magnitude,
            g_junction_net.to("MPa").magnitude,
            places=10,
        )

    def test_weight_fractions_from_universe_with_multibead_b1_and_b0(self):
        """MMT weight fractions from a Universe must use chain functionality, not bead degree."""
        from pylimer_tools_cpp import Universe

        # Construct a network with A4 (type 2), multi-bead B2 (type 1),
        # multi-bead B1 (type 3), and multi-bead B0 solvent (type 4).
        # 2 crosslinkers (f=4): 8 A sites; 7 bonds formed => p = 7/8 = 0.875
        # 3 B2 chains (3 beads each): 6 B sites
        # 2 B1 chains (3 beads each): 2 B sites
        # 2 B0 chains (3 beads each): 0 B sites
        # r = 8 / 8 = 1.0, b2 = 6 / 8 = 0.75
        net = Universe(10, 10, 10)
        atom_ids = (
            [1, 2] + list(range(10, 19)) + list(range(20, 26)) + list(range(30, 36))
        )
        atom_types = [2, 2] + [1] * 9 + [3] * 6 + [4] * 6
        n_atoms = len(atom_ids)
        net.add_atoms(
            atom_ids,
            atom_types,
            [float(i) for i in range(n_atoms)],
            [0.0] * n_atoms,
            [0.0] * n_atoms,
            [0] * n_atoms,
            [0] * n_atoms,
            [0] * n_atoms,
        )
        # Internal chain bonds (all chains have 3 beads => internal bead degree == 2)
        bonds_from = [10, 11, 13, 14, 16, 17, 20, 21, 23, 24, 30, 31, 33, 34]
        bonds_to = [11, 12, 14, 15, 17, 18, 21, 22, 24, 25, 31, 32, 34, 35]
        # Crosslinker bonds:
        # xlink 1 gets 4 bonds (f=4): B2 ends 10, 13, 16 and B1 end 20
        # xlink 2 gets 3 bonds: B2 ends 12, 15 and B1 end 23 (B2 end 18 left unreacted)
        bonds_from += [1, 12, 1, 15, 1, 1, 2]
        bonds_to += [10, 2, 13, 2, 16, 20, 23]
        net.add_bonds(bonds_from, bonds_to)
        net.set_masses({1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0})

        # Confirm bead degree is 2 for B1 (type 3) and B0 (type 4)
        bead_func = net.determine_functionality_per_type()
        self.assertEqual(bead_func[3], 2)
        self.assertEqual(bead_func[4], 2)

        # Compute MMT weight fractions directly from the Universe (functionality_per_type=None)
        w_sol_auto = compute_weight_fraction_of_soluble_material(
            net, crosslinker_type=2
        )
        w_dang_auto = compute_weight_fraction_of_dangling_chains(
            net, crosslinker_type=2
        )
        w_bb_auto = compute_weight_fraction_of_backbone(net, crosslinker_type=2)

        # Mass conservation must hold exactly
        self.assertAlmostEqual(w_sol_auto + w_dang_auto + w_bb_auto, 1.0, places=10)

        # Compare against explicit analytical parameters
        wf = net.compute_weight_fractions()
        explicit_params = {
            "crosslinker_type": 2,
            "functionality_per_type": {1: 2, 2: 4, 3: 1, 4: 0},
            "weight_fractions": wf,
            "r": 1.0,
            "p": 7.0 / 8.0,
            "b2": 0.75,
        }
        w_sol_ref = compute_weight_fraction_of_soluble_material(None, **explicit_params)
        w_dang_ref = compute_weight_fraction_of_dangling_chains(None, **explicit_params)
        w_bb_ref = compute_weight_fraction_of_backbone(None, **explicit_params)

        self.assertAlmostEqual(w_sol_auto, w_sol_ref, places=10)
        self.assertAlmostEqual(w_dang_auto, w_dang_ref, places=10)
        self.assertAlmostEqual(w_bb_auto, w_bb_ref, places=10)

        # Verify solvent contributes its full weight fraction to w_sol and 0 to w_dang & w_bb
        _, beta = compute_miller_macosko_probabilities(r=1.0, p=7.0 / 8.0, f=4, b2=0.75)
        self.assertGreater(w_sol_auto, wf[4] + wf[3] * beta)
        # If B1/B0 were mistakenly treated as B2 (beta^2 instead of beta and 1.0),
        # w_sol would be strictly smaller since beta < 1.
        wrong_w_sol = compute_weight_fraction_of_soluble_material(
            None,
            crosslinker_type=2,
            functionality_per_type={1: 2, 2: 4, 3: 2, 4: 2},
            weight_fractions=wf,
            r=1.0,
            p=7.0 / 8.0,
            b2=0.75,
        )
        self.assertGreater(w_sol_auto, wrong_w_sol)

        # Partial functionality_per_type={2: 4} with network must infer B2/B1/B0
        # without KeyError and without counting solvent (type 4) in nu
        self.assertAlmostEqual(
            compute_weight_fraction_of_soluble_material(
                net, crosslinker_type=2, functionality_per_type={2: 4}
            ),
            w_sol_ref,
            places=10,
        )
        self.assertAlmostEqual(
            compute_weight_fraction_of_dangling_chains(
                net, crosslinker_type=2, functionality_per_type={2: 4}
            ),
            w_dang_ref,
            places=10,
        )
        self.assertAlmostEqual(
            compute_weight_fraction_of_backbone(
                net, crosslinker_type=2, functionality_per_type={2: 4}
            ),
            w_bb_ref,
            places=10,
        )
        unit_style = UnitStyleFactory().get_unit_style("si")
        g_auto_partial, _, g_anm_partial, _ = compute_modulus_decomposition(
            network=net,
            unit_style=unit_style,
            crosslinker_type=2,
            functionality_per_type={2: 4},
        )
        g_auto_full, _, g_anm_full, _ = compute_modulus_decomposition(
            network=net,
            unit_style=unit_style,
            crosslinker_type=2,
            functionality_per_type={1: 2, 2: 4, 3: 1, 4: 0},
        )
        self.assertAlmostEqual(
            g_auto_partial.to("MPa").magnitude,
            g_auto_full.to("MPa").magnitude,
            places=10,
        )
        self.assertAlmostEqual(
            g_anm_partial.to("MPa").magnitude,
            g_anm_full.to("MPa").magnitude,
            places=10,
        )

    def test_pure_monofunctional_b2_zero_and_ulp_p_max(self) -> None:
        """Pure A_f + B_1 (b2 == 0) and 1-ULP saturated r > 1 networks must not crash."""
        from pylimer_tools_cpp import Universe

        unit_style = UnitStyleFactory().get_unit_style("si")
        ureg = unit_style.get_underlying_unit_registry()

        alpha, beta = compute_miller_macosko_probabilities(r=1.0, p=0.8, f=4, b2=0.0)
        self.assertEqual((alpha, beta), (1.0, 1.0))
        self.assertTrue(np.isinf(predict_gelation_point(r=1.0, f=4, b2=0.0)))

        with self.assertRaises(ValueError):
            compute_miller_macosko_probabilities(r=1.0, p=0.8, f=4, b2=-0.1)

        with self.assertRaises(ValueError):
            predict_p_from_w_sol(
                w_sol=0.5,
                network=None,
                crosslinker_type=2,
                functionality_per_type={2: 4, 3: 1},
                weight_fractions={2: 0.2, 3: 0.8},
                r=1.0,
                b2=0.0,
            )

        # Build a pure A4 + B1 Universe (1 crosslinker, 4 2-bead B1 chains)
        b1_net = Universe(10, 10, 10)
        b1_ids = [1, 10, 11, 12, 13, 14, 15, 16, 17]
        b1_types = [2] + [3] * 8
        b1_net.add_atoms(
            b1_ids,
            b1_types,
            [float(i) for i in range(len(b1_ids))],
            [0.0] * len(b1_ids),
            [0.0] * len(b1_ids),
            [0] * len(b1_ids),
            [0] * len(b1_ids),
            [0] * len(b1_ids),
        )
        b1_net.add_bonds([10, 12, 14, 16, 1, 1, 1, 1], [11, 13, 15, 17, 10, 12, 14, 16])
        b1_net.set_masses({2: 1.0, 3: 1.0})

        self.assertAlmostEqual(
            compute_weight_fraction_of_soluble_material(b1_net, crosslinker_type=2),
            1.0,
            places=10,
        )
        self.assertAlmostEqual(
            compute_weight_fraction_of_dangling_chains(b1_net, crosslinker_type=2),
            0.0,
            places=10,
        )
        self.assertAlmostEqual(
            compute_weight_fraction_of_backbone(b1_net, crosslinker_type=2),
            0.0,
            places=10,
        )
        g_ph, g_ent, _, _ = compute_modulus_decomposition(
            network=b1_net, unit_style=unit_style, crosslinker_type=2
        )
        self.assertAlmostEqual(g_ph.to("MPa").magnitude, 0.0, places=10)
        self.assertAlmostEqual(g_ent.to("MPa").magnitude, 0.0, places=10)

        # Saturated r > 1 network where (10 / 7) / 3 exceeds 1 / 2.1 by 1 ULP in IEEE-754
        r_ulp = 21.0 / 10.0
        p_ulp = (10.0 / 7.0) / 3.0
        alpha_u, beta_u = compute_miller_macosko_probabilities(r=r_ulp, p=p_ulp, f=3)
        self.assertTrue(0.0 <= alpha_u <= 1.0)
        self.assertTrue(0.0 <= beta_u <= 1.0)

        # Explicit f=4 override on a network whose max observed crosslinker degree is 3
        partial_deg_net = Universe(10, 10, 10)
        partial_deg_net.add_atoms(
            [1, 2, 10, 11, 12, 13],
            [2, 2, 1, 1, 1, 1],
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            [0.0] * 6,
            [0.0] * 6,
            [0] * 6,
            [0] * 6,
            [0] * 6,
        )
        # xlink 1 has 3 bonds (10, 11, 12), xlink 2 has 3 bonds (10, 11, 13) -> max deg = 3
        partial_deg_net.add_bonds([1, 10, 1, 11, 1, 2], [10, 2, 11, 2, 12, 13])
        partial_deg_net.set_masses({1: 1.0, 2: 1.0})
        g_ent_f4 = compute_entanglement_modulus(
            0.5 * ureg("MPa"), network=partial_deg_net, crosslinker_type=2, f=4
        )
        g_ent_ref = compute_entanglement_modulus(
            0.5 * ureg("MPa"), p=6.0 / 8.0, r=8.0 / 8.0, f=4, b2=1.0
        )
        self.assertAlmostEqual(
            g_ent_f4.to("MPa").magnitude, g_ent_ref.to("MPa").magnitude, places=10
        )


if __name__ == "__main__":
    unittest.main()
