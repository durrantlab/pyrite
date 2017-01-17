import unittest
import os
import sys

import numpy as np
import scipy
import scoria


class InformationTests(unittest.TestCase):
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

    def test_get_filename(self):
        """
        Tests the getting of filenames.
        """
        expected_filename = ["./scoria/sample_files/3_mol_test.pdb"]
        self.assertEqual(self.mol.get_filename(), expected_filename)

    def test_get_remarks(self):
        """
        Tests the getting of remarks.
        """
        expected_remarks = [" This is a test file."]
        self.assertEqual(self.mol.get_remarks(), expected_remarks)

    def test_get_center_of_mass(self):
        """
        Tests the determination of the center of mass.
        """
        expected_center = [-21.57417476, 52.60448475, -17.17966907]
        center_of_mass = self.mol.get_center_of_mass()
        self.assertAlmostEqual(center_of_mass[0], expected_center[0], self.accuracy)
        self.assertAlmostEqual(center_of_mass[1], expected_center[1], self.accuracy)
        self.assertAlmostEqual(center_of_mass[2], expected_center[2], self.accuracy)

    @unittest.skip("Needs final value")
    def test_get_atom_information(self):
        """
        Tests the atom information.
        """
        expected_atom_inf = [] # WRite final value here
        atom_inf = self.mol.get_atom_information()
        self.assertListEqual(atom_inf, expected_atom_inf)


    def test_get_coordinates(self):
        """
        Tests that the coordinates are returned as expected.
        """
        coordinates = self.mol.get_coordinates()
        self.assertEqual(len(coordinates), 12)

        atom_zero = coordinates[0]
        expected_atom = [10.000, 10.000, 10.000]
        self.assertAlmostEqual(atom_zero[0], expected_atom[0], self.accuracy)
        self.assertAlmostEqual(atom_zero[1], expected_atom[1], self.accuracy)
        self.assertAlmostEqual(atom_zero[2], expected_atom[2], self.accuracy)

    def test_get_bonds(self):
        """
        Tests that the bond values are returned.
        """
        bonds = self.mol.get_bonds()

        # Checking that it is the correct depth, length and value
        self.assertEqual(len(bonds), 12)
        self.assertEqual(len(bonds[0]), 12)
        self.assertEqual(bonds[0][0], 0)
        self.assertEqual(bonds[2][3], 1)

    def test_get_geometric_center(self):
        """
        Tests that the geometric center is able to be calculated properly.
        """
        expected_center = [-21.75291634, 52.4852562, -17.28250122]
        geo_center = self.mol.get_geometric_center()

        self.assertAlmostEqual(geo_center[0], expected_center[0], self.accuracy)
        self.assertAlmostEqual(geo_center[1], expected_center[1], self.accuracy)
        self.assertAlmostEqual(geo_center[2], expected_center[2], self.accuracy)

    def test_get_total_number_of_atoms(self):
        """
        Tests that the number of atoms returned is correct.
        """
        number_of_atoms = self.mol.get_total_number_of_atoms()
        self.assertEqual(number_of_atoms, 12)

    def test_get_total_mass(self):
        """
        Tests that the total mass is returned correctly.
        """
        expected_mass = 156.104
        total_mass = self.mol.get_total_mass()
        self.assertAlmostEqual(total_mass, expected_mass, self.accuracy)

    # Depreciated? And needs skip for dependencies
    @unittest.skip("Deprecated Function?")
    def test_get_heirarchy(self):
        """
        Tests that the hierarchy can be set.
        """
        hierarchy = self.mol.get_hierarchy()
        # Assertation here

    ## Testing Setters

    def test_set_filename(self):
        """
        Tests the setting of filenames.
        """
        self.mol.set_filename("OtherFile.txt")
        self.assertEqual(self.mol.get_filename(), ["OtherFile.txt"])

    def test_set_remarks(self):
        """
        Tests the setting of remarks.
        """
        set_remarks = ["TEST REMARK"]
        self.mol.set_remarks(set_remarks)
        self.assertEqual(self.mol.get_remarks(), set_remarks)

    @unittest.skip("Needs test written")
    def test_set_coordinates(self):
        """
        Tests that the coordinates can be set
        """
        coordinates = []
        self.mol.set_coordinates(coordinates)
        # Assertation here

    @unittest.skip("Needs test written")
    def test_set_atomic_information(self):
        """
        Tests that the atom information can be set.
        """
        atom_inf = {}
        self.mol.set_atom_information(atom_inf)
        # Assertation here

    # Add skip for dependencies
    @unittest.skip("Needs test written")
    def test_set_bonds(self):
        """
        Tests that the atom information can be set.
        """
        bonds = {}
        self.mol.set_bonds(bonds)
        # Assertation here

    @unittest.skip("Needs test written")
    def test_set_coordinate_undo_point(self):
        """
        Tests that the coordinate undo point can be set.
        """
        coord_undo = {}
        self.mol.set_coordinates_undo_point(coord_undo)
        # Assertation here

    # Depreciated? And needs skip for dependencies
    @unittest.skip("Needs test written")
    def test_set_heirarchy(self):
        """
        Tests that the hierarchy can be set.
        """
        hierarchy = {}
        self.mol.set_hierarchy(hierarchy)
        # Assertation here


    ## Testing Functions

    # The bounding box, having several parameters, should have some
    # comprehensive tests written for it.
    @unittest.skip("Needs test written")
    def test_get_default_bounding_box(self):
        """
        Tests that the bounding box can be calculated.
        """
        bounding_box = self.mol.get_bounding_box(None, 0.0, 0)
        # Assertation here

    # Similar to the bounding box tests, we need to check that all
    # parameters work properly.
    @unittest.skip("Needs test written")
    def test_get_bounding_sphere(self):
        """
        Tests that the bounding sphere can be calculated.
        """
        bounding_sphere = self.mol.get_bounding_sphere(None, 0.0, 0)
        # Assertation here

    @unittest.skip("Needs test written")
    def test_get_constants(self):
        """
        Tests that the constants returned are as expected
        """
        constants = self.mol.get_constants(self)
        # Assertation here

    # For the 'belongs' tests, we need one index of each category
    # And each should be tested for redundancy
    @unittest.skip("Needs correct test indexes")
    def test_belongs_to_protein(self):
        """
        Tests that indices are proteins
        """
        self.assertTrue(self.mol.belongs_to_protein(0))
        self.assertFalse(self.mol.belongs_to_protein(1))
        self.assertFalse(self.mol.belongs_to_protein(2))

    @unittest.skip("Needs correct test indexes")
    def test_belongs_to_dna(self):
        """
        Tests that indices are DNA
        """
        self.assertTrue(self.mol.belongs_to_dna(0))
        self.assertFalse(self.mol.belongs_to_dna(1))
        self.assertFalse(self.mol.belongs_to_dna(2))

    @unittest.skip("Needs correct test indexes")
    def test_belongs_to_RNA(self):
        """
        Tests that indices are RNA
        """
        self.assertTrue(self.mol.belongs_to_rna(0))
        self.assertFalse(self.mol.belongs_to_rna(1))
        self.assertFalse(self.mol.belongs_to_rna(2))

    @unittest.skip("Needs test written")
    def test_assign_elements_from_atom_names(self):
        """
        Tests the assignment of elements from the atom names.
        """
        # Assertion here, pre assignment
        self.mol.assign_elements_from_atom_names()
        # Assertion here, post assignment

    @unittest.skip("Needs test written")
    def test_assign_masses(self):
        """
        Tests the assignment of masses.
        """
        # Assertion here, pre assignment
        self.mol.assign_masses()
        # Assertion here, post assignment

    @unittest.skip("Needs test written")
    def test_serial_reindex(self):
        """
        Tests the reindexing of the serial field.
        """
        # Assertion here, pre assignment
        self.mol.serial_reindex()
        # Assertion here, post assignment

    @unittest.skip("Needs test written")
    def test_resseq_reindex(self):
        """
        Tests the reindexing of the resseq field.
        """
        # Assertion here, pre assignment
        self.mol.resseq_reindex()
        # Assertion here, post assignment

    @unittest.skip("Needs test written")
    def test_define_molecule_chain_residue_spherical_boundaries(self):
        """
        Tests the reindexing of the serial field.
        """
        # Assertion here, pre assignment
        self.mol.define_molecule_chain_residue_spherical_boundaries()
        # Assertion here, post assignment


if __name__ == '__main__':
    unittest.main()