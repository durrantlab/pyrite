"""
Copyright (c) 2017 Jacob Durrant. MIT license. Please see LICENSE.txt for full details.
"""
from __future__ import absolute_import
import warnings
import unittest
from . import InformationTests as IT
from . import FileIOTests as FIOT
from . import GeometryTests as GT
from . import ManipulationTests as MT
from . import OtherMoleculeTests as OMT
from . import SelectionTests as ST

warnings.filterwarnings("ignore", category=DeprecationWarning)

class UnitTests(object):
    """
    Unit testing object for scoria.
    """
    def __init__(self):
        """
        Initalizes the unit tests.
        """
        self._suite = unittest.TestSuite()
        self._runner = unittest.TextTestRunner()

    # Running Suite

    def run(self):
        """
        Runs the currently queued suite of tests.
        """
        self._runner.run(self._suite)

    def run_all(self):
        """
        Quickly runs all unit tests.
        """
        self.add_all_tests()
        self.run()

    # Add specific module tests

    def add_all_tests(self):
        """
        Adds all available tests to the suite.
        """
        self.add_information_tests()
        self.add_fileio_tests()
        self.add_geometry_tests()
        self.add_selection_tests()
        self.add_manipulation_tests()
        self.add_other_molecule_tests()

    def add_information_tests(self):
        """
        Adds the information tests.
        """
        information_tests = unittest.makeSuite(IT.InformationTests)
        self._suite.addTests(information_tests)

    def add_fileio_tests(self):
        """
        Adds the information tests.
        """
        fileio_tests = unittest.makeSuite(FIOT.FileIOTests)
        self._suite.addTests(fileio_tests)

    def add_geometry_tests(self):
        """
        Adds the information tests.
        """
        tests = unittest.makeSuite(GT.GeometryTests)
        self._suite.addTests(tests)
        
    def add_manipulation_tests(self):
        """
        Adds the information tests.
        """
        tests = unittest.makeSuite(MT.ManipulationTests)
        self._suite.addTests(tests)

    def add_other_molecule_tests(self):
        """
        Adds the information tests.
        """
        tests = unittest.makeSuite(OMT.OtherMoleculeTests)
        self._suite.addTests(tests)

    def add_selection_tests(self):
        """
        Adds the information tests.
        """
        tests = unittest.makeSuite(ST.SelectionsTests)
        self._suite.addTests(tests)
        
        