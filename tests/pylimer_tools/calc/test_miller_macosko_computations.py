import copy
import itertools
import os
import sys
import unittest

from pylimer_tools.calc.miller_macosko_theory import (
    compute_miller_macosko_probabilities,
    compute_modulus_decomposition,
    compute_weight_fraction_of_soluble_material,
    compute_weight_fractions_and_probabilities,
    predict_gelation_point,
    predict_number_density_of_junction_points,
    compute_weight_fraction_of_backbone,
    compute_weight_fraction_of_dangling_chains,
    predict_number_density_of_network_strands,
    predict_shear_modulus,
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
            compute_stoichiometric_imbalance(
                self.testUniverse, 2, effective=True),
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
            compute_stoichiometric_imbalance(
                self.testUniverse, 2, ignore_types=[3]),
            ((3 * 3) / ((5 * 2) / (5 / 3))),
        )

    def test_extent_of_reaction(self):
        self.assertAlmostEqual(
            1.0, compute_extent_of_reaction(
                self.emptyUniverse, 2))
        self.assertAlmostEqual(
            5.0 / 6.0, compute_extent_of_reaction(self.testUniverse, 2)
        )

    def test_gelation_point_prediction(self):
        self.assertAlmostEqual(1, predict_gelation_point(1, 2))
        self.assertAlmostEqual(1, predict_gelation_point(1, 2, 2))

    def test_shear_modulus_prediction(self):
        self.assertRaises(
            ValueError,
            lambda: predict_shear_modulus(
                network=self.emptyUniverse, crosslinker_type=2, unit_style=None
            ),
        )
        self.assertRaises(
            ValueError,
            lambda: predict_shear_modulus(
                network=self.emptyUniverse, unit_style=None),
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

    def test_predict_number_density_of_junction_points(self):
        self.testUniverse.set_masses({1: 1, 2: 1})
        self.assertAlmostEqual(
            1.5, compute_stoichiometric_imbalance(self.testUniverse, 2)
        )
        # TODO: find literature motiviation for results fo the functions
        self.assertAlmostEqual(
            predict_number_density_of_junction_points(
                self.testUniverse, 2), 0.0
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
            predict_number_density_of_network_strands(
                self.testUniverse, 2), 0.0
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
        self.assertDictEqual(
            weight_fractions, {
                1: 1 - 3.0 / 8.0, 2: 3.0 / 8.0})
        test_universe_copy = copy.copy(self.testUniverse)
        test_universe_copy.remove_atoms([1, 2, 3, 4, 5, 6])
        self.assertTrue(test_universe_copy.get_nr_of_atoms() == 2)
        self.assertDictEqual(
            compute_weight_fractions(test_universe_copy), {
                1: 1.0 / 2.0, 2: 1.0 / 2.0}
        )

        params = {
            "crosslinker_type": 2,
            "functionality_per_type": {1: 2, 2: 3},
            "r": 1.0,
            "p": 0.95,
            "weight_fractions": {1: 0.85, 2: 0.15},
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

    def test_soluble_weight_fraction_measurement(self):
        self.testUniverse.set_masses({1: 1, 2: 1})
        self.assertEqual(
            measure_weight_fraction_of_soluble_material(
                self.emptyUniverse), 0.0
        )
        self.assertEqual(
            measure_weight_fraction_of_soluble_material(
                self.testUniverse, rel_tol=0),
            0.0,
        )
        self.assertEqual(
            measure_weight_fraction_of_soluble_material(
                self.testUniverse, abs_tol=10000
            ),
            1.0,
        )
        self.assertEqual(
            measure_weight_fraction_of_soluble_material(
                self.testUniverse), 1 / 8
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
            NotImplementedError,
            lambda: compute_weight_fraction_of_soluble_material(
                self.testUniverse, 2, functionality_per_type={1: 1, 2: 3}
            ),
        )
        self.assertRaises(
            ValueError,
            lambda: compute_weight_fraction_of_soluble_material(
                self.testUniverse, 7),
        )
        self.saturatedTestUniverse.set_masses({1: 1, 2: 1})

        res_tuple = compute_weight_fractions_and_probabilities(
            self.saturatedTestUniverse, 2
        )
        expected_tuple = ({1: 0.85, 2: 0.15},
                          0.111111111111111, 0.111111111111111)
        self.assertEqual(len(res_tuple), len(expected_tuple))
        for i in range(len(res_tuple)):
            if isinstance(res_tuple[i], int) or isinstance(
                    res_tuple[i], float):
                self.assertAlmostEqual(res_tuple[i], expected_tuple[i])
            elif isinstance(res_tuple[i], dict):
                for key in res_tuple[i].keys():
                    self.assertAlmostEqual(
                        res_tuple[i][key], expected_tuple[i][key])
            else:
                raise ValueError(
                    "Expected integer, float or dict for comparison")

        self.assertAlmostEqual(
            0.010699588477366243,
            compute_weight_fraction_of_soluble_material(
                self.saturatedTestUniverse, 2),
        )

    # def testProbabilityCalculations(self):
    #     self.assertRaises(
    #         ValueError, lambda: compute_miller_macosko_probabilities(0.9, 2, 2))
    #     self.assertRaises(
    # ValueError, lambda: compute_miller_macosko_probabilities(0.1, 0.9, 3))

    def test_backbone_weight_fraction_calculations(self):
        self.assertEqual(
            0, compute_weight_fraction_of_backbone(
                self.emptyUniverse, 2))
        self.assertEqual(
            0.0, compute_weight_fraction_of_dangling_chains(
                self.emptyUniverse, 2)
        )
        self.assertEqual(
            1, compute_weight_fraction_of_soluble_material(
                self.emptyUniverse, 2)
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
            compute_weight_fraction_of_dangling_chains(
                self.saturatedTestUniverse, 2),
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
            nu=4.69218e25 *
            (unit_style.get_underlying_unit_registry()("meter") ** -3),
            temperature=298 * unit_style.get_underlying_unit_registry()("kelvin"),
        )

        alpha, beta = compute_miller_macosko_probabilities(
            r=1.0, p=0.95, f=4.0)
        self.assertAlmostEqual(alpha, 0.0983588, places=5)

        self.assertAlmostEqual(g_anm.to("MPa").magnitude, 0.1930520, places=5)
        self.assertAlmostEqual(g_pnm.to("MPa").magnitude, 0.0965260, places=5)
        self.assertAlmostEqual(
            g_mmt_entanglement.to("MPa").magnitude, 0.190576, places=5
        )
        self.assertAlmostEqual(g_mmt_phantom.to(
            "MPa").magnitude, 0.0777123, places=5)

        # these in turn require further investigation into the involvement of r
        g_mmt_phantom, g_mmt_entanglement, g_anm, g_pnm = compute_modulus_decomposition(
            network=None,
            unit_style=unit_style,
            crosslinker_type=2,
            r=1.3,
            p=0.6465,
            f=4,
            nu=1.25981e25 *
            (unit_style.get_underlying_unit_registry()("meter") ** -3),
            temperature=298 * unit_style.get_underlying_unit_registry()("kelvin"),
        )

        self.assertAlmostEqual(g_anm.to("MPa").magnitude, 0.051833, places=5)
        self.assertAlmostEqual(g_pnm.to("MPa").magnitude, 0.025916, places=5)
        self.assertAlmostEqual(
            g_mmt_entanglement.to("MPa").magnitude, 0.05801087, places=5
        )
        self.assertAlmostEqual(
            g_mmt_phantom.to("MPa").magnitude,
            0.00492674,
            places=5)

    def test_probability_calculations(self):
        # most already in other tests
        alpha, beta = compute_miller_macosko_probabilities(r=1.0, p=0.85, f=5)
        self.assertAlmostEqual(0.282074, alpha, places=5)
        alpha, beta = compute_miller_macosko_probabilities(r=1.0, p=0.85, f=6)
        self.assertAlmostEqual(0.278715, alpha, places=5)
        self.assertRaises(
            ValueError, lambda: compute_miller_macosko_probabilities(
                r=1.0, p=2, f=7)
        )


if __name__ == "__main__":
    unittest.main()
