import os
import sys
import unittest

import numpy as np

from pylimer_tools.calc.mehp_utilities import (
    compute_cycle_rank, compute_effective_nr_density_of_junctions,
    compute_effective_nr_density_of_network, compute_mean_universe_volume,
    compute_topological_factor, predict_shear_modulus)
from pylimer_tools.calc.structure_analysis import (
    compute_crossLinker_conversion,
    compute_effective_crossLinker_functionalities,
    compute_effective_crossLinker_functionality,
    compute_mean_end_to_end_distances, compute_mean_end_to_end_vectors,
    measure_weight_fraction_of_backbone,
    measure_weight_fraction_of_dangling_chains)
from pylimer_tools_cpp.pylimer_tools_cpp import MoleculeType, Universe

if __name__ == '__main__':
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMEHPAnalysisFunctions(UniverseUsingTestCase):

    def test_weight_fraction_calculations(self):
        self.assertEqual(
            (0.0, 0.0), measure_weight_fraction_of_dangling_chains(self.emptyUniverse, 2))
        # empty weight -> empty weight fraction
        self.testUniverse.setMasses({1: 0, 2: 0})
        self.assertEqual(
            (0.0, 0.25), measure_weight_fraction_of_backbone(self.testUniverse, 2))
        self.assertEqual(
            1.0, measure_weight_fraction_of_backbone(self.testUniverse, 2))
        # non-empty weights
        self.testUniverse.setMasses({1: 1, 2: 0})
        self.assertTrue(self.testUniverse.getNrOfAtoms() > 0)
        self.assertEqual(self.testUniverse.getMasses(), {1: 1, 2: 0})
        all_chains = self.testUniverse.getChainsWithCrosslinker(2)
        self.assertEqual(all_chains[2].getType(), MoleculeType.DANGLING_CHAIN)
        self.assertEqual(
            (0.2, 0.25), measure_weight_fraction_of_dangling_chains(self.testUniverse, crossLinker_type=2))

    def test_crossLinker_functionality_calculation(self):
        self.assertCountEqual(
            [], compute_effective_crossLinker_functionalities(self.emptyUniverse, 2))
        self.assertSequenceEqual(
            [0, 2, 3], compute_effective_crossLinker_functionalities(self.testUniverse, 2))
        self.assertEqual(
            5.0 / 3.0, compute_effective_crossLinker_functionality(self.testUniverse, 2))
        self.assertEqual(
            5.0 / 3.0 / 3.0, compute_crossLinker_conversion(self.testUniverse, 2, 3))

    def test_mean_end_to_end_computation(self):
        self.assertCountEqual([], compute_mean_end_to_end_vectors([], 2))
        self.assertDictEqual({
            '1-2-3-6-7': 1.0,
            '5-6-7': 1.0
        }, compute_mean_end_to_end_distances([self.testUniverse], 2))

    def test_mean_universe_volume(self):
        self.assertRaises(NotImplementedError, lambda: compute_mean_universe_volume(
            [self.testUniverse, self.testUniverseSmall]))

    def test_effective_nr_density_of_junction_calculation(self):
        self.assertIsNone(compute_effective_nr_density_of_junctions([]))
        # Border cases
        self.assertEqual(
            0.0, compute_effective_nr_density_of_junctions([self.testUniverse], 0, 0, crossLinker_type=None))
        self.assertEqual(
            0.0, compute_effective_nr_density_of_junctions([self.testUniverse], 1000, crossLinker_type=2))
        self.assertEqual(
            0.0, compute_effective_nr_density_of_junctions([self.emptyUniverse], 1000, crossLinker_type=2))
        # Other border
        # 3 junctions, volume of 1
        self.assertEqual(
            3.0 / self.testUniverse.getVolume(),
            compute_effective_nr_density_of_junctions([self.testUniverse], 0,
                                                      crossLinker_type=2, min_num_effective_strands=0))
        # actual calc: 6 & 7 are active, 4 not
        self.assertEqual(
            2.0 / self.testUniverse.getVolume(),
            compute_effective_nr_density_of_junctions([self.testUniverse], abs_tol=None, rel_tol=0,
                                                      crossLinker_type=2, min_num_effective_strands=2))
        self.assertEqual(
            2.0 / self.testUniverse.getVolume(),
            compute_effective_nr_density_of_junctions([self.testUniverse], 0,
                                                      crossLinker_type=2, min_num_effective_strands=2))

    def test_effective_nr_density_of_network_calculation(self):
        self.assertIsNone(compute_effective_nr_density_of_network([]))
        self.assertEqual(3, len(self.testUniverse.getMolecules(2)))
        # Border cases
        self.assertEqual(0.0, compute_effective_nr_density_of_network(
            [self.testUniverse], None, 10, crossLinker_type=2))
        self.assertEqual(
            0.0, compute_effective_nr_density_of_network([self.testUniverse], 100, 100, crossLinker_type=2))
        self.assertEqual(
            0.0, compute_effective_nr_density_of_network([self.testUniverse], 1000, 1, crossLinker_type=2))
        # actual calc: we got 2 active strands in a Volume of 1
        self.assertEqual(
            2.0 / self.testUniverse.getVolume(),
            compute_effective_nr_density_of_network([self.testUniverse], 0, 2,
                                                    crossLinker_type=2))

    def test_cycle_rank_calculation(self):
        self.assertEqual(1, compute_cycle_rank(None, 1, 0))
        self.assertEqual(0, compute_cycle_rank(None, 1, 1))
        self.assertEqual(-1, compute_cycle_rank(None, 0, 1))
        universe = Universe(10, 10, 10)
        universe = self.addAtomBondData(
            universe, self.testAtoms, self.testBonds)
        # test basic exception thrown when specifying the wrong arguments
        self.assertRaises(ValueError, lambda: compute_cycle_rank(
            networks=[universe], crossLinker_type=None))
        self.assertRaises(
            ValueError, lambda: compute_cycle_rank(networks=[universe], nu=1, crossLinker_type=None))
        # same nr of active strands as junctions
        self.assertEqual(
            0.0, compute_cycle_rank(networks=[universe], nu=None, mu=None, abs_tol=1, rel_tol=1, crossLinker_type=2))
        # other system
        self.assertEqual(
            1 / (10 * 10 * 10), compute_cycle_rank(networks=[self.saturatedTestUniverse], nu=None, mu=None, abs_tol=1, rel_tol=1, crossLinker_type=2))

    def test_topological_factor_computation(self):
        self.assertEqual(
            1 + 1.0 / 3.0, compute_topological_factor([self.testUniverse], 2, b=1))
        bond_lengths = []
        for m in self.testUniverse.getMolecules(2):
            bond_lengths.extend(m.computeBondLengths())
        self.assertEqual(1, np.mean(bond_lengths))
        self.assertEqual(
            0.5485762961986437, compute_topological_factor([self.testUniverse], 2))
        # larger system
        # g = self.saturatedTestUniverse.getUnderlyingGraph()
        # igraph.plot(g, vertex_label=g.vs["name"], vertex_color=["green" if n["type"] == 2 else "red" for n in g.vs],
        # target="large_test.png", vertex_label_dist=1)
        self.assertEqual(0.7249043914053506, compute_topological_factor(
            [self.saturatedTestUniverse], 2))

    def test_shear_modulus_prediction(self):
        self.assertEqual(0.0, predict_shear_modulus(
            [self.emptyUniverse], crossLinker_type=2))
        self.assertEqual(0.003624521957026753, predict_shear_modulus(
            [self.saturatedTestUniverse], crossLinker_type=2))


if __name__ == '__main__':
    unittest.main()
