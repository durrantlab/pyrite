import unittest
import os
import sys

import numpy as np
import scipy
import scoria


class GeometryTests(unittest.TestCase):
    """
    Base Test Suite
    """

    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        self.mol = scoria.Molecule()
        self.accuracy = 4

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        self.mol = None

        # Remove all files from the tmp folder?

    ### Tests

    def test_get_angle_between_three_points(self):
        """
        Empty test.
        """
        one = np.array([1.0, 0.0, 0.0])
        two = np.array([0.0, 0.0, 0.0])
        three = np.array([0.0, 0.0, 1.0])

        radians_right = self.mol.get_angle_between_three_points(one, two, three)
        degrees_right = np.degrees(radians_right)

        radians_less_one = self.mol.get_angle_between_three_points(two, three, one)
        degrees_less_one = np.degrees(radians_less_one)

        radians_less_two = self.mol.get_angle_between_three_points(three, one, two)
        degrees_less_two = np.degrees(radians_less_two)

        self.assertAlmostEqual(degrees_less_one, degrees_less_two, self.accuracy)
        self.assertAlmostEqual(degrees_less_one, 45.0, self.accuracy)
        self.assertAlmostEqual(degrees_right, 90.0, self.accuracy)

    def test_get_dihedral_angle(self):
        """
        Empty test.
        """
        one = np.array([0.0, -1.0, 1.0])
        two = np.array([0.0, -1.0, 0.0])
        three = np.array([0.0, 1.0, 0.0])
        four = np.array([1.0, 1.0, 0.0])

        radians = self.mol.get_dihedral_angle(one, two, three, four)
        degrees = np.degrees(radians)

        self.assertAlmostEqual(degrees, 90.0, self.accuracy)

    def test_is_planar(self):
        """
        Empty test.
        """
        one = np.array([0.0, -1.0, 1.0])
        two = np.array([0.0, -1.0, 0.0])
        three = np.array([0.0, 1.0, 0.0])
        four = np.array([1.0, 1.0, 0.0])

        self.assertFalse(self.mol.is_planar(one, two, three, four))

        four = np.array([0.0, 1.0, 1.0])

        self.assertTrue(self.mol.is_planar(one, two, three, four))

    # We need to figure out how this functions
    def test_get_planarity_deviation(self):
        """
        Empty test.
        """
        one = np.array([0.0, 100.0, 100.0])
        two = np.array([0.0, 100.0, 0.0])
        three = np.array([0.0, 0.0, 100.0])
        four = np.array([5.0, 0.0, 0.0])

        deviation = self.mol.get_planarity_deviation(one, two, three, four)

        self.assertAlmostEqual(deviation, 5.0, self.accuracy)
