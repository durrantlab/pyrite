import unittest
import os
import sys

import numpy as np
import scipy
import scoria


class FileIOTests(unittest.TestCase):
    """
    Base Test Suite
    """

    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        if not os.path.exists("./scoria_tests_tmp"):
            os.mkdir("./scoria_tests_tmp")

        self.mol = scoria.Molecule("PDB", "./scoria/sample_files/3_mol_test.pdb")
        self.accuracy = 4

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.mol = None

        # Remove all files from the tmp folder?

    ### Tests
    # Testing Getters

    def test_nothing(self):
        """
        Empty test.
        """
        pass
