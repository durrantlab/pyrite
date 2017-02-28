from __future__ import absolute_import
import unittest
import os
import sys

import numpy as np
import scipy
import scoria


class OtherMoleculeTests(unittest.TestCase):
    """
    Base Test Suite
    """

    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        self.mol = scoria.Molecule()

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.mol = None

        # Remove all files from the tmp folder?

    ### Tests

    @unittest.skip("Needs test written")
    def test_realignment(self):
        """
        Empty test.
        """
        
