#!/usr/bin/env python
import unittest

from pylimer_tools_cpp import Universe


class TestEntitiesNew(unittest.TestCase):

    def test_universe(self):
        universe = Universe(10, 10, 10)
        self.assertIsInstance(universe, Universe)


if __name__ == '__main__':
    unittest.main()
