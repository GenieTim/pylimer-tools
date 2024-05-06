import random
import unittest

import numpy as np
from numpy import linalg

from pylimer_tools_cpp import Box


class TestMEHPAssumptions(unittest.TestCase):
    def test_box_offset_distribution(self):
        """
        This test tests, whether the assumptions made to determine
        the distribution of the box offsets when introducing a slip-link
        is correct.
        """
        box = Box(10, 10, 10)
        coords1 = np.array([5, 2.5, 0.])
        coords2 = np.array([25, 7.5, 0.])
        offset = box.get_offset(coords2 - coords1)
        self.assertAlmostEqual(offset[0], -20.)
        self.assertAlmostEqual(offset[1], 0.)
        self.assertAlmostEqual(offset[2], 0.)

        for _ in range(10):
            # attempt 10 times
            # find random coords in the box
            slip_link_coords = np.array(
                [random.random() * 10, random.random() * 10, random.random() * 10])
            #
            dist_1 = slip_link_coords - coords1
            box.apply_pbc(dist_1)
            dist_2 = coords2 - slip_link_coords
            box.apply_pbc(dist_2)

            # we have two choices, and want to
            reference_norm = min(linalg.norm(dist_1 + offset) + linalg.norm(
                dist_2), linalg.norm(dist_1) + linalg.norm(dist_2 + offset))

            # with this choice done, this test here
            # should now check that no other combination gets smaller


if __name__ == '__main__':
    unittest.main()
