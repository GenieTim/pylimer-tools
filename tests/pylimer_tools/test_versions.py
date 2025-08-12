#!/usr/bin/env python

import importlib.metadata as importlib_metadata
import unittest

import pylimer_tools as pt
import pylimer_tools_cpp as ptc


class TestVersions(unittest.TestCase):

    def test_pylimer_tools_version(self):
        self.assertIsNotNone(pt.__version__)

    def test_pylimer_tools_cpp_version(self):
        self.assertIsNotNone(ptc.__version__)

    def test_equivalent_versions(self):
        self.assertEqual(pt.__version__, ptc.__version__)
        self.assertEqual(importlib_metadata.version("pylimer-tools"), pt.__version__)


if __name__ == "__main__":
    unittest.main()
