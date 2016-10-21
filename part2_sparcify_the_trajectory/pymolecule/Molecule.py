import numpy
from FileIO import FileIO
from AtomsAndBonds import AtomsAndBonds
from Selections import Selections
from Manipulation import Manipulation
from Information import Information
from OtherMolecules import OtherMolecules
from Geometry import Geometry


class Molecule: # here's the actual Molecule class
    """Loads, saves, and manupulates molecuar models. The main pymolecule
    class."""

    def __init__ (self):
        """Initializes the variables of the Molecule class."""

        self.fileio = FileIO(self)
        self.atoms_and_bonds = AtomsAndBonds(self)
        self.selections = Selections(self)
        self.manipulation = Manipulation(self)
        self.information = Information(self)
        self.other_molecule = OtherMolecules(self)
        self.geometry = Geometry(self)

    # Information methods
    ### Aliases ###
    # Gets
    def get_coordinates(self):
        return self.information.get_coordinates()

    def get_filename(self):
        return self.information.get_filename()

    def get_remarks(self):
        return self.information.get_remarks()

    def get_atom_information(self):
        return self.information.get_atom_information()

    def get_coordinates_undo_point(self):
        return self.information.get_coordinates_undo_point()

    def get_bonds(self):
        return self.information.get_bonds()

    def get_hierarchy(self):
        return self.information.get_hierarchy()

    def get_constants(self):
        return self.information.get_constants()

    def get_center_of_mass(self, selection = None):
        return self.information.get_center_of_mass(selection)

    def get_geometric_center(self, selection = None):
        return self.information.get_geometric_center(selection)

    def get_total_mass(self, selection = None):
        return self.information.get_total_mass(selection)

    def get_total_number_of_atoms(self, selection = None):
        return self.information.get_total_number_of_atoms(selection)

    def get_total_number_of_heavy_atoms(self, selection = None):
        return self.information.get_total_number_of_heavy_atoms(selection)

    def get_bounding_box(self, selection = None, padding = 0.0):
        return self.information.get_bounding_box(selection, padding)

    def get_bounding_sphere(self, selection = None, padding = 0.0):
        return self.information.get_bounding_sphere(selection, padding)

    # Set
    def set_filename(self, filename):
        self.information.set_filename(filename)

    def set_remarks(self, remarks):
        self.information.set_remarks(remarks)

    def set_atom_information(self, atom_information):
        self.information.set_atom_information(atom_information)

    def set_coordinates(self, coordinates):
        self.information.set_coordinates(coordinates)

    def set_coordinates_undo_point(self, coordinates_undo_point):
        self.information.set_coordinates_undo_point(coordinates_undo_point)

    def set_bonds(self, bonds):
        self.information.set_bonds(bonds)

    def set_hierarchy(self, hierarchy):
        self.information.set_hierarchy(hierarchy)

    # Information functions
    def assign_masses(self):
        self.information.assign_masses()

    def assign_elements_from_atom_names(self, selection = None):
        self.information.assign_elements_from_atom_names(selection)

    def define_molecule_chain_residue_spherical_boundaries(self):
        self.information.define_molecule_chain_residue_spherical_boundaries()

    def serial_reindex(self):
        self.information.serial_reindex()

    def resseq_reindex(self):
        self.information.resseq_reindex()

    # File I/O class methods
    def load_pym_into(self, filename):
        self.fileio.load_pym_into(filename)

    def load_pdb_into(self, filename, bonds_by_distance = True,
                      serial_reindex = True, resseq_reindex = False):

        self.fileio.load_pdb_into(
            filename, bonds_by_distance, serial_reindex, resseq_reindex
        )

    def load_pdb_into_using_file_object(self, file_obj,
                                        bonds_by_distance = True,
                                        serial_reindex = True,
                                        resseq_reindex = False):

        self.fileio.load_pdb_into_using_file_object(
            file_obj, bonds_by_distance, serial_reindex, resseq_reindex
        )

    def load_pdbqt_into(self, filename, bonds_by_distance = False,
                      serial_reindex = True, resseq_reindex = False):

        self.fileio.load_pdbqt_into(
            filename, bonds_by_distance, serial_reindex, resseq_reindex
        )

    def load_pdbqt_into_using_file_object(self, file_obj,
                                        bonds_by_distance = False,
                                        serial_reindex = True,
                                        resseq_reindex = False):

        self.fileio.load_pdbqt_into_using_file_object(
            file_obj, bonds_by_distance, serial_reindex, resseq_reindex
        )

    def save_pym(self, filename, save_bonds = False, save_filename = False,
                 save_remarks = False, save_hierarchy = False,
                 save_coordinates_undo_point = False):

        self.fileio.save_pym(
            filename, save_bonds, save_filename, save_remarks,
            save_hierarchy, save_coordinates_undo_point
        )

    def save_pdb(self, filename = "", serial_reindex = True,
                 resseq_reindex = False, return_text = False):

        self.fileio.save_pdb(
            filename, serial_reindex, resseq_reindex, return_text
        )

    # Atoms and Bonds class methods
    def get_number_of_bond_partners_of_element(self, atom_index, the_element):

        return self.atoms_and_bonds.get_number_of_bond_partners_of_element(
            atom_index, the_element
        )

    def get_index_of_first_bond_partner_of_element(self, atom_index,
                                                   the_element):

        return self.atoms_and_bonds.get_index_of_first_bond_partner_of_element(
            atom_index, the_element
        )

    def create_bonds_by_distance(self, remove_old_bond_data = True,
                                 delete_excessive_bonds = True):
        self.atoms_and_bonds.create_bonds_by_distance(
            remove_old_bond_data, delete_excessive_bonds
        )

    def delete_bond(self, index1, index2):
        self.atoms_and_bonds.delete_bond(index1, index2)

    def add_bond(self, index1, index2, order = 1):
        self.atoms_and_bonds.add_bond(index1, index2, order)

    def delete_atom(self, index):
        self.atoms_and_bonds.delete_atom(index)

    def add_atom(self, record_name = "ATOM", serial = 1, name = "X",
                 resname = "XXX", chainid = "X", resseq = 1, occupancy = 0.0,
                 tempfactor = 0.0, charge = '', element = "X",
                 coordinates = numpy.array([0.0, 0.0, 0.0]), autoindex = True):

        self.atoms_and_bonds.add_atom(
            record_name, serial, name, resname, chainid, resseq, occupancy,
            tempfactor, charge, element, coordinates, autoindex
        )

    # Selections class
    def get_molecule_from_selection(self, selection, serial_reindex = True,
                                    resseq_reindex = False):

        return self.selections.get_molecule_from_selection(
            selection, serial_reindex, resseq_reindex
        )

    def select_atoms(self, selection_criteria):
        return self.selections.select_atoms(selection_criteria)

    def select_atoms_in_bounding_box(self, bounding_box):
        return self.selections.select_atoms_in_bounding_box(bounding_box)

    def select_branch(self, root_atom_index, directionality_atom_index):
        return self.selections.select_branch(
            root_atom_index, directionality_atom_index
        )

    def select_all_atoms_bound_to_selection(self, selections):
        return self.selections.select_all_atoms_bound_to_selection(selections)

    def select_atoms_from_same_molecule(self, selection):
        return self.selections.select_atoms_from_same_molecule(root_atom_index)

    def selections_of_constituent_molecules(self):
        return self.selections.selections_of_constituent_molecules()

    def select_atoms_near_other_selection(self, selection, cutoff):
        return self.selections.select_atoms_near_other_selection(
            selection, cutoff
        )

    def select_atoms_in_same_residue(self, selection):
        return self.selections.select_atoms_in_same_residue(selection)

    def invert_selection(self, selection):
        return self.selections.invert_selection(selection)

    def select_all(self):
        return self.selections.select_all()

    def select_close_atoms_from_different_molecules(self, other_mol, cutoff,
                                                    pairwise_comparison = True,
                                                    terminate_early = False):
        return self.selections.select_close_atoms_from_different_molecules(
            other_mol, cutoff, pairwise_comparison, terminate_early
        )

    def selections_of_chains(self):
        return self.selections.selections_of_chains()

    def selections_of_residues(self):
        return self.selections.selections_of_residues()

    # Manipulation class
    def set_atom_location(self, atom_index, new_location):
        return self.manipulation.set_atom_location(atom_index, new_location)

    def coordinate_undo(self):
        self.manipulation.coordinate_undo()

    def translate_molecule(self, delta):
        self.manipulation.translate_molecule(delta)

    def rotate_molecule_around_a_line_between_points(self, line_point1,
                                                     line_point2, rotate):

        self.manipulation.rotate_molecule_around_a_line_between_points(
            line_point1, line_point2, rotate
        )

    def rotate_molecule_around_a_line_between_atoms(self, line_point1_index,
                                                    line_point2_index, rotate):

        self.manipulation.rotate_molecule_around_a_line_between_atoms(
            line_point1_index, line_point2_index, rotate
        )

    def rotate_molecule_around_pivot_point(self, pivot, thetax,
                                           thetay, thetaz):

        self.manipulation.rotate_molecule_around_pivot_point(
            pivot, thetax, thetay, thetaz
        )

    def rotate_molecule_around_pivot_atom(self, pivot_index, thetax,
                                          thetay, thetaz):

        self.manipulation.rotate_molecule_around_pivot_atom(
            pivot_index, thetax, thetay, thetaz
        )

    # Geometry class
    def get_angle_between_three_points(self, pt1, pt2, pt3):
        return self.geometry.get_angle_between_three_points(pt1, pt2, pt3)

    def get_dihedral_angle(self, pt1, pt2, pt3, pt4):
        return self.geometry.get_dihedral_angle(pt1, pt2, pt3, pt4)

    def get_planarity_deviation(self, pt1, pt2, pt3, pt4):
        return self.geometry.get_planarity_deviation(pt1, pt2, pt3, pt4)

    def is_planar(self, pt1, pt2, pt3, pt4, planarity_cutoff = 0.2):
        return self.geometry.is_planar(pt1, pt2, pt3, pt4, planarity_cutoff)

    # Other molecule class
    def get_other_molecule_aligned_to_this(self, other_mol, tethers):
        # Add Weight Matrix
        return self.other_molecule.get_other_molecule_aligned_to_this(
            other_mol, tethers
        )

    def get_distance_to_another_molecule(self, other_molecule,
                                         pairwise_comparison = True):

        return self.other_molecule.get_distance_to_another_molecule(
            other_moelcule, pairwise_comparison
        )

    def get_rmsd_equivalent_atoms_specified(self, other_mol, tethers):

        return self.other_molecule.get_rmsd_equivalent_atoms_specified(
            other_mol, tethers
        )

    def get_rmsd_order_dependent(self, other_mol):
        return self.other_molecule.get_rmsd_order_dependent(other_mol)

    def get_rmsd_heuristic(self, other_mol):
        return self.other_molecule.get_rmsd_heuristic(other_mol)

    def steric_clash_with_another_molecule(self, other_mol, cutoff,
                                           pairwise_comparison = True):

        return self.other_molecule.steric_clash_with_another_molecule(
            other_mol, cutoff, pairwise_comparison
        )

    def merge_with_another_molecule(self, other_molecule):
        return self.other_molecule.merge_with_another_molecule(other_molecule)

    ######## Supporting functions ########

    def numpy_structured_array_remove_field(self, narray, field_names):
        """Removes a specific field name from a structured numpy array.

            Args::
                narray -- A structured numpy array.
                field_names -- A list of strings, where each string is one of
                    the field names of narray.

            Returns:
                A structured numpy array identical to narray, but with the
                    field names in field_names removed.

        """

        # surprised this doesn't come with numpy

        # now remove the coordinates from the atom_information object to save
        # memory
        names = list(narray.dtype.names)
        for f in field_names: names.remove(f)
        return narray[names]

    def __is_number(self, s):
        """Determines whether or not a string represents a number.

            Args::
                s -- A string (e.g., "5.4").

            Returns:
                A boolean, whether or not the string can be represented by a
                    float.

        """

        try:
            float(s)
            return True
        except ValueError:
            return False

    def copy(self):
        """Returns an exact copy (pymolecule.Molecule) of this Molecule object.
        Undo points are NOT copied.

            Returns:
                A pymolecule.Molecule, containing to the same atomic
                    information as this pymolecule.Molecule object.

        """

        new_molecule = Molecule()
        new_molecule.set_filename(self.get_filename()[:])
        new_molecule.set_remarks(self.get_remarks()[:])
        new_molecule.set_atom_information(self.get_atom_information().copy())
        new_molecule.set_coordinates(self.get_coordinates().copy())

        if not self.get_bonds() is None:
            new_molecule.set_bonds(self.get_bonds().copy())
        else:
            new_molecule.set_bonds(None)

        new_molecule.set_hierarchy(copy.deepcopy(self.get_hierarchy()))

        return new_molecule
