from Molecule import Molecule
import cStringIO as StringIO
import os
from pymolecule import dumbpy as numpy


class Test:
    """A class for testing all pymolecule functions."""

    mol = None

    def test_all(self):
        """Test all pymolecule functions."""

        print "Creating directory to store temporary files."

        if not os.path.exists("./pymolecule_tests_tmp"):
            os.mkdir("./pymolecule_tests_tmp")

        self.test_file_io()
        #self.test_information()
        #self.test_selection()
        #self.test_manipulation()
        #self.test_other_molecules()
        #self.test_atoms_and_bonds()
        #self.test_geometry()

    def test_file_io(self):
        """Test the functions in FileIO."""

        file_io_filename = "./pymolecule_tests_tmp/file_io_test"

        print "FileIO Functions"
        print "    load_pdb_into_using_file_object()"
        self.mol = Molecule()
        self.mol.load_pdb_into_using_file_object(
            StringIO.StringIO(self.pdb_file_str()),
            True,
            True,
            True
        )

        print "    save_pdb()"
        self.mol.save_pdb(file_io_filename + ".pdb", True, True, False)

        print "    save_pym()"
        self.mol.save_pym(file_io_filename + ".pym", True, True, True, True, True)

        print "    load_pdbqt_into_using_file_object()"
        self.mol = Molecule()
        self.mol.load_pdb_into_using_file_object(
            StringIO.StringIO(self.pdbqt_file_str()),
            False,  # no bonds by distance, because pdbqt file has unrecognized atom types.
            True,
            True
        )

        open(file_io_filename + ".pdbqt", 'w').write(self.pdbqt_file_str())

        print "    load_pdbqt_into()"
        self.mol = Molecule()
        self.mol.load_pdbqt_into(file_io_filename + ".pdbqt", False, True, True)

        print "    load_pym_into()"
        self.mol = Molecule()
        self.mol.load_pym_into(file_io_filename + ".pym")
    
        print "    load_pdb_into()"
        self.mol = Molecule()
        self.mol.load_pdb_into(file_io_filename + ".pdb", True, True, True)

    def test_information(self):
        """Test the functions in Information"""

        print "Information Functions"
        print "    get_filename()"
        print "        Filename: " + self.mol.get_filename()

        print "    get_remarks()"
        print "        " + str(self.mol.get_remarks())

        print "    assign_elements_from_atom_names()"
        self.mol.assign_elements_from_atom_names()

        print "    assign_masses()"
        self.mol.assign_masses()

        print "    get_center_of_mass()"
        print "        " + str(self.mol.get_center_of_mass())

        print "    get_atom_information()"
        print "        " + str(self.mol.get_atom_information()[0])

        print "    get_coordinates()"
        print "        " + str(self.mol.get_coordinates()[0])

        print "    get_bonds()"
        print "        " + str(self.mol.get_bonds()[0][:10])

        print "    get_geometric_center()"
        print "        " + str(self.mol.get_geometric_center())

        print "    get_total_number_of_atoms()"
        print "        " + str(self.mol.get_total_number_of_atoms())

        print "    get_total_number_of_heavy_atoms()"
        print "        " + str(self.mol.get_total_number_of_heavy_atoms())

        print "    get_total_mass()"
        print "        " + str(self.mol.get_total_mass())

        print "    get_bounding_box()"
        print "        " + str(self.mol.get_bounding_box())

        print "    get_bounding_sphere()"
        print "        " + str(self.mol.get_bounding_sphere())

        print "    get_constants()"
        print "        " + str(self.mol.get_constants())

        print "    belongs_to_protein()"
        print "        " + str(self.mol.belongs_to_protein(1))

        print "    belongs_to_dna()"
        print "        " + str(self.mol.belongs_to_dna(1))

        print "    belongs_to_rna()"
        print "        " + str(self.mol.belongs_to_rna(1))

        print "    set_filename()"
        self.mol.set_filename("test_name.pdb")
        print "        " + str(self.mol.get_filename())

        print "    set_remarks()"
        self.mol.set_remarks("Test remark.")
        print "        " + str(self.mol.get_remarks())

        print "    set_coordinates()"
        coors = self.mol.get_coordinates()
        coors[:,0] = 999.999
        self.mol.set_coordinates(coors)
        print "        " + str(self.mol.get_coordinates()[0])

        print "    set_atom_information()"
        inf = self.mol.get_atom_information()
        inf["resname"] = " TST "
        self.mol.set_atom_information(inf)
        print "        " + str(self.mol.get_atom_information()[0])

        print "    set_bonds()"
        bnds = self.mol.get_bonds()
        bnds[:,:10] = 1
        self.mol.set_bonds(bnds)
        print "        " + str(self.mol.get_bonds()[0][:10])

        print "    serial_reindex()"
        self.mol.serial_reindex()
        
        print "    resseq_reindex()"
        self.mol.resseq_reindex()
        
        print "    set_coordinates_undo_point()"
        self.mol.set_coordinates_undo_point(
            self.mol.get_coordinates() + 10
        )
        
        print "    define_molecule_chain_residue_spherical_boundaries()"
        self.mol.define_molecule_chain_residue_spherical_boundaries()
        
        print "    get_hierarchy()"
        print "        " + str(self.mol.get_hierarchy()['residues']['indices']['VAL-32-A'])

        print "    set_hierarchy()"
        self.mol.set_hierarchy(
            self.mol.get_hierarchy()
        )

    def test_selection(self):
        """Test the functions in Selections."""

        print "Selection Functions"

        # Get new molecule... clean slate.
        self.mol = Molecule()
        self.mol.load_pdb_into_using_file_object(
            StringIO.StringIO(self.pdb_file_str()),
            True,
            True,
            True
        )

        print "    select_all()"
        print "        Atoms in selection: " + str(
            len(
                self.mol.select_all()
            )
        )

        print "    select_atoms()"
        sel = self.mol.select_atoms({"resname_stripped": "TRP"})
        print "        Atoms in selection: " + str(len(sel))

        print "    invert_selection()"
        print "        Atoms in selection: " + str(len(self.mol.invert_selection(sel)))
        
        print "    select_all_atoms_bound_to_selection()"
        print "        Atoms in selection: " + str(len(self.mol.select_all_atoms_bound_to_selection(sel)))
        
        print "    select_atoms_near_other_selection()"
        print "        Atoms in selection: " + str(len(self.mol.select_atoms_near_other_selection(sel, 8.0)))
        
        print "    select_atoms_in_same_residue()"
        print "        Atoms in selection: " + str(len(self.mol.select_atoms_in_same_residue([1])))
        
        print "    select_atoms_from_same_molecule()"
        print "        Atoms in selection: " + str(len(self.mol.select_atoms_from_same_molecule([1])))
        
        print "    select_atoms_from_same_molecule()"
        print "        Atoms in selection: " + str(len(self.mol.select_atoms_from_same_molecule([1])))
        

        print "    select_atoms_in_bounding_box()"
        print "        Atoms in selection: " + str(
            len(
                self.mol.select_atoms_in_bounding_box(
                    self.mol.get_bounding_box(range(15), 5.0)
                )
            )
        )

        print "    select_branch()"  # Get a side chain.
        print "        Atoms in selection: " + str(len(self.mol.select_branch(1, 4)))

        print "    selections_of_constituent_molecules()" 
        print "        Atoms in selection: " + str(self.mol.selections_of_constituent_molecules())

        print "    selections_of_chains()" 
        print "        Chains mapped to indices: " + str(self.mol.selections_of_chains().keys())

        print "    selections_of_residues()" 
        print "        Residues mapped to indices: " + str(self.mol.selections_of_residues().keys()[:5])

        print "    select_close_atoms_from_different_molecules()"
        sels = self.mol.select_close_atoms_from_different_molecules(self.mol, 1.0)
        print "        Atoms in mol1 selection: " + str(len(sels[0]))
        print "        Atoms in mol2 selection: " + str(len(sels[1]))

        print "    get_molecule_from_selection()"
        sel = self.mol.select_branch(1, 4)
        print "        Atoms in selection: " + str(len(sel))
        mol2 = self.mol.get_molecule_from_selection(sel)
        mol2.save_pdb("./pymolecule_tests_tmp/save_selection.pdb", True, True, False)

    def test_manipulation(self):
        """Test the functions in Manipulation."""

        manip_filename = "./pymolecule_tests_tmp/manipulation_test"

        print "Manipulation Functions"
        
        print "    set_coordinates_undo_point()"
        self.mol.set_coordinates_undo_point(
            self.mol.get_coordinates()
        )
        
        print "    translate_molecule()"
        self.mol.translate_molecule(numpy.array([10.0, 10.0, 10.0]))
        self.mol.save_pdb(manip_filename + "1.pdb", False, False, False)
        
        print "    set_atom_location()"
        self.mol.set_atom_location(1, numpy.array([10.0, 10.0, 10.0]))
        self.mol.save_pdb(manip_filename + "2.pdb", False, False, False)
        
        print "    rotate_molecule_around_a_line_between_points()"
        self.mol.rotate_molecule_around_a_line_between_points(
            numpy.array([0.0, 0.0, 0.0]),
            numpy.array([10.0, 0.0, 0.0]),
            45.0
        )
        self.mol.save_pdb(manip_filename + "3.pdb", False, False, False)
        
        print "    rotate_molecule_around_a_line_between_atoms()"
        self.mol.rotate_molecule_around_a_line_between_atoms(1, 2, 45.0)
        self.mol.save_pdb(manip_filename + "4.pdb", False, False, False)
        
        print "    rotate_molecule_around_pivot_point()"
        self.mol.rotate_molecule_around_pivot_point(numpy.array([0.0, 0.0, 0.0]), 45, 45, 45)
        self.mol.save_pdb(manip_filename + "5.pdb", False, False, False)
        
        print "    rotate_molecule_around_pivot_atom()"
        self.mol.rotate_molecule_around_pivot_atom(5, 45, 45, 45)
        self.mol.save_pdb(manip_filename + "6.pdb", False, False, False)
        
        print "    coordinate_undo()"
        self.mol.coordinate_undo()
        self.mol.save_pdb(manip_filename + "7.pdb", False, False, False)
        
    def test_other_molecules(self):
        """Test the functions in OtherMolecules"""

        print "OtherMolecules Functions"

        # Make two molecules to play with.
        mol1 = Molecule()
        mol1.load_pdb_into_using_file_object(
            StringIO.StringIO(self.pdb_file_str()),
            True,
            True,
            True
        )

        mol2 = Molecule()
        mol2.load_pdb_into_using_file_object(
            StringIO.StringIO(self.pdb_file_str()),
            True,
            True,
            True
        )

        mol2.translate_molecule(numpy.array([10.0, 10.0, 10.0]))

        print "    steric_clash_with_another_molecule()"
        print "        " + str(mol1.steric_clash_with_another_molecule(mol2, 5.0, False))
        print "        " + str(mol1.steric_clash_with_another_molecule(mol2, 5.0, True))

        print "    get_distance_to_another_molecule()"
        print "        " + str(mol1.get_distance_to_another_molecule(mol2, False))
        print "        " + str(mol1.get_distance_to_another_molecule(mol2, True))

        print "    get_rmsd_order_dependent()"
        print "        " + str(mol1.get_rmsd_order_dependent(mol2))

        print "    get_rmsd_heuristic()"
        print "        " + str(mol1.get_rmsd_heuristic(mol2))

        print "    get_rmsd_equivalent_atoms_specified()"
        tethers = [numpy.array([1, 1]), numpy.array([2, 2]), numpy.array([3, 3])]
        print "        " + str(mol1.get_rmsd_equivalent_atoms_specified(mol2, tethers))

        print "    merge_with_another_molecule()"
        merged_mol = mol1.merge_with_another_molecule(mol2)
        merged_mol.save_pdb("./pymolecule_tests_tmp/merged.pdb", False, False, False)

        print "    get_other_molecule_aligned_to_this()"
        aligned_mol = mol1.get_other_molecule_aligned_to_this(mol2, tethers)
        print "        New RMSD: " + str(mol1.get_rmsd_order_dependent(aligned_mol))

    def test_atoms_and_bonds(self):
        """Test the functions in AtomsAndBonds."""

        self.mol = Molecule()
        
        print "AtomsAndBonds Functions"

        print "    add_atom()"
        self.mol.add_atom()
        self.mol.add_atom()
        self.mol.add_atom()
        self.mol.add_atom()
        self.mol.add_atom()

        self.mol.set_coordinates(numpy.array(
            [[0.0, 0.0, 0.0],
             [1.5, 0.0, 0.0],
             [1.5, 1.5, 0.0],
             [3.0, 1.5, 0.0],
             [999.99, 2.0, 0.0]]
        ))

        print "    delete_atom()"
        self.mol.delete_atom(4)

        print "    add_bond()"
        self.mol.add_bond(0, 1, 2)

        print "    delete_bond()"
        self.mol.delete_bond(0, 1)

        print "    create_bonds_by_distance()"
        self.mol.create_bonds_by_distance()

        print "    get_number_of_bond_partners_of_element()"
        print "        " + str(self.mol.get_number_of_bond_partners_of_element(0, "X"))

        print "    get_index_of_first_bond_partner_of_element()"
        print "        " + str(self.mol.get_index_of_first_bond_partner_of_element(1, "X"))

    def test_geometry(self):
        """Test the functions in Geometry."""

        print "Geometry Functions"

        coors = self.mol.get_coordinates()

        print "    get_angle_between_three_points()"
        print "        " + str(self.mol.get_angle_between_three_points(coors[0], coors[1], coors[2]))

        print "    get_dihedral_angle()"
        print "        " + str(self.mol.get_dihedral_angle(coors[0], coors[1], coors[2], coors[3]))

        print "    is_planar()"
        print "        " + str(self.mol.is_planar(coors[0], coors[1], coors[2], coors[3]))

        print "    get_planarity_deviation()"
        print "        " + str(self.mol.get_planarity_deviation(coors[0], coors[1], coors[2], coors[3]))

    def pdb_file_str(self):
        """Get a PDB-formatted string.
        
            Returns:
                A string, PDB formatted.
        """

        return """REMARK This is a remark.
ATOM      1  N   GLN A  52      42.237  16.800  35.823  1.00 12.04           N  
ATOM      2  CA  GLN A  52      41.667  17.015  34.477  1.00 10.21           C  
ATOM      3  C   GLN A  52      40.148  17.010  34.446  1.00 10.33           C  
ATOM      4  O   GLN A  52      39.564  17.006  33.382  1.00  9.88           O  
ATOM      5  CB  GLN A  52      42.228  15.967  33.522  1.00  9.36           C  
ATOM      6  CG  GLN A  52      43.734  16.095  33.352  1.00  9.40           C  
ATOM      7  CD  GLN A  52      44.170  17.346  32.622  1.00  8.47           C  
ATOM      8  OE1 GLN A  52      45.312  17.791  32.785  1.00 12.47           O  
ATOM      9  NE2 GLN A  52      43.296  17.915  31.849  1.00  8.79           N  
ATOM     10  N   SER A  53      39.488  17.025  35.604  1.00 11.01           N  
ATOM     11  CA  SER A  53      38.026  16.949  35.611  1.00 12.75           C  
ATOM     12  C   SER A  53      37.340  18.091  34.862  1.00 12.38           C  
ATOM     13  O   SER A  53      36.209  17.937  34.426  1.00 14.69           O  
ATOM     14  CB  SER A  53      37.462  16.890  37.026  1.00 13.59           C  
ATOM     15  OG  SER A  53      37.854  18.006  37.781  1.00 17.90           O  
ATOM     16  N   ASP A  54      38.012  19.228  34.730  1.00 11.65           N  
ATOM     17  CA  ASP A  54      37.455  20.341  33.968  1.00 11.85           C  
ATOM     18  C   ASP A  54      37.742  20.271  32.474  1.00 10.23           C  
ATOM     19  O   ASP A  54      37.168  21.045  31.722  1.00 11.10           O  
ATOM     20  CB  ASP A  54      37.947  21.667  34.531  1.00 13.04           C  
ATOM     21  CG  ASP A  54      37.464  21.907  35.937  1.00 16.08           C  
ATOM     22  OD1 ASP A  54      36.280  21.638  36.211  1.00 19.99           O  
ATOM     23  OD2 ASP A  54      38.207  22.360  36.811  1.00 21.76           O  
ATOM     24  N   PHE A  55      38.617  19.367  32.044  1.00  9.13           N  
ATOM     25  CA  PHE A  55      38.972  19.219  30.642  1.00  8.35           C  
ATOM     26  C   PHE A  55      37.915  18.485  29.848  1.00  8.87           C  
ATOM     27  O   PHE A  55      37.392  17.475  30.303  1.00 10.67           O  
ATOM     28  CB  PHE A  55      40.300  18.470  30.543  1.00  8.52           C  
ATOM     29  CG  PHE A  55      40.730  18.127  29.138  1.00  7.26           C  
ATOM     30  CD1 PHE A  55      41.050  19.110  28.229  1.00  7.47           C  
ATOM     31  CD2 PHE A  55      40.847  16.803  28.753  1.00  8.45           C  
ATOM     32  CE1 PHE A  55      41.470  18.763  26.957  1.00  8.40           C  
ATOM     33  CE2 PHE A  55      41.266  16.455  27.495  1.00  8.90           C  
ATOM     34  CZ  PHE A  55      41.584  17.426  26.589  1.00  9.37           C  
ATOM     35  N   SER A  56      37.615  18.984  28.653  1.00  8.47           N  
ATOM     36  CA  SER A  56      36.702  18.341  27.727  1.00  8.75           C  
ATOM     37  C   SER A  56      37.409  18.135  26.390  1.00  8.22           C  
ATOM     38  O   SER A  56      37.676  19.116  25.681  1.00  8.56           O  
ATOM     39  CB  SER A  56      35.476  19.231  27.540  1.00  9.47           C  
ATOM     40  OG  SER A  56      34.454  18.554  26.839  1.00 11.56           O  
ATOM     41  N   PRO A  57      37.729  16.888  26.033  1.00  8.33           N  
ATOM     42  CA  PRO A  57      38.409  16.649  24.754  1.00  8.76           C  
ATOM     43  C   PRO A  57      37.622  17.191  23.585  1.00  8.53           C  
ATOM     44  O   PRO A  57      36.400  17.105  23.551  1.00 10.15           O  
ATOM     45  CB  PRO A  57      38.525  15.124  24.677  1.00  9.69           C  
ATOM     46  CG  PRO A  57      38.403  14.640  26.059  1.00 11.00           C  
ATOM     47  CD  PRO A  57      37.508  15.628  26.772  1.00  9.30           C  
ATOM     48  N   TYR A  58      38.339  17.747  22.619  1.00  7.67           N  
ATOM     49  CA  TYR A  58      37.722  18.174  21.383  1.00  8.58           C  
ATOM     50  C   TYR A  58      37.720  16.916  20.539  1.00 10.08           C  
ATOM     51  O   TYR A  58      37.954  15.802  21.051  1.00 12.68           O  
ATOM     52  CB  TYR A  58      38.459  19.389  20.783  1.00  7.74           C  
ATOM     53  CG  TYR A  58      37.649  20.058  19.693  1.00  7.14           C  
ATOM     54  CD1 TYR A  58      36.490  20.749  19.989  1.00  7.63           C  
ATOM     55  CD2 TYR A  58      38.011  19.951  18.354  1.00  7.35           C  
ATOM     56  CE1 TYR A  58      35.725  21.321  18.986  1.00  6.97           C  
ATOM     57  CE2 TYR A  58      37.237  20.520  17.348  1.00  7.25           C  
ATOM     58  CZ  TYR A  58      36.095  21.188  17.671  1.00  6.13           C  
ATOM     59  OH  TYR A  58      35.360  21.745  16.658  1.00  7.63           O  
ATOM     60  N   ILE A  59      37.333  17.049  19.295  1.00 12.35           N  
ATOM     61  CA  ILE A  59      37.115  15.901  18.448  1.00 12.48           C  
ATOM     62  C   ILE A  59      38.089  15.820  17.302  1.00 11.76           C  
ATOM     63  O   ILE A  59      38.750  16.803  16.917  1.00 13.39           O  
ATOM     64  CB  ILE A  59      35.694  15.926  17.889  1.00 13.52           C  
ATOM     65  CG1 ILE A  59      35.445  17.222  17.112  1.00 15.44           C  
ATOM     66  CG2 ILE A  59      34.704  15.759  19.016  1.00 15.64           C  
ATOM     67  CD1 ILE A  59      34.256  17.154  16.212  1.00 17.64           C  
ATOM     68  N   GLU A  60      38.180  14.611  16.777  1.00 10.67           N  
ATOM     69  CA  GLU A  60      38.661  14.409  15.428  1.00 10.55           C  
ATOM     70  C   GLU A  60      37.481  14.459  14.452  1.00  8.69           C  
ATOM     71  O   GLU A  60      36.332  14.278  14.839  1.00 10.07           O  
ATOM     72  CB  GLU A  60      39.495  13.129  15.350  1.00 12.80           C  
ATOM     73  CG  GLU A  60      40.798  13.316  16.137  1.00 17.70           C  
ATOM     74  CD  GLU A  60      41.613  12.057  16.350  1.00 21.64           C  
ATOM     75  OE1 GLU A  60      41.570  11.155  15.484  1.00 23.57           O  
ATOM     76  OE2 GLU A  60      42.308  11.976  17.395  1.00 25.13           O  
ATOM     77  N   ILE A  61      37.793  14.733  13.206  1.00  6.90           N  
ATOM     78  CA  ILE A  61      36.779  14.950  12.200  1.00  7.13           C  
ATOM     79  C   ILE A  61      37.084  14.056  11.030  1.00  7.33           C  
ATOM     80  O   ILE A  61      38.219  13.998  10.572  1.00  8.79           O  
ATOM     81  CB  ILE A  61      36.757  16.443  11.793  1.00  6.83           C  
ATOM     82  CG1 ILE A  61      36.419  17.333  12.993  1.00  6.74           C  
ATOM     83  CG2 ILE A  61      35.749  16.663  10.695  1.00  6.97           C  
ATOM     84  CD1 ILE A  61      36.548  18.813  12.744  1.00  7.60           C  
ATOM     85  N   ASP A  62      36.066  13.372  10.530  1.00  6.96           N  
ATOM     86  CA  ASP A  62      36.254  12.425   9.448  1.00  7.34           C  
ATOM     87  C   ASP A  62      36.040  13.011   8.060  1.00  6.96           C  
ATOM     88  O   ASP A  62      35.232  13.916   7.872  1.00  6.52           O  
ATOM     89  CB  ASP A  62      35.294  11.244   9.583  1.00  8.39           C  
ATOM     90  CG  ASP A  62      35.411  10.521  10.881  1.00 11.97           C  
ATOM     91  OD1 ASP A  62      36.506  10.494  11.474  1.00 15.65           O  
ATOM     92  OD2 ASP A  62      34.418   9.952  11.370  1.00 17.75           O  
ATOM     93  N   LEU A  63      36.758  12.454   7.084  1.00  6.71           N  
ATOM     94  CA  LEU A  63      36.396  12.571   5.682  1.00  6.81           C  
ATOM     95  C   LEU A  63      35.017  11.951   5.504  1.00  7.10           C  
ATOM     96  O   LEU A  63      34.647  11.040   6.252  1.00  7.91           O  
ATOM     97  CB  LEU A  63      37.396  11.824   4.813  1.00  7.22           C  
ATOM     98  CG  LEU A  63      38.828  12.357   4.856  1.00  7.73           C  
ATOM     99  CD1 LEU A  63      39.771  11.385   4.148  1.00  8.36           C  
ATOM    100  CD2 LEU A  63      38.933  13.737   4.253  1.00  8.58           C  
ATOM    101  N   PRO A  64      34.242  12.400   4.526  1.00  7.30           N  
ATOM    102  CA  PRO A  64      32.886  11.862   4.349  1.00  7.38           C  
ATOM    103  C   PRO A  64      32.877  10.416   3.868  1.00  7.74           C  
ATOM    104  O   PRO A  64      33.404  10.094   2.797  1.00  8.70           O  
ATOM    105  CB  PRO A  64      32.270  12.793   3.299  1.00  8.12           C  
ATOM    106  CG  PRO A  64      33.462  13.342   2.535  1.00  7.13           C  
ATOM    107  CD  PRO A  64      34.540  13.481   3.572  1.00  7.55           C  
ATOM    108  N   SER A  65      32.257   9.546   4.658  1.00  8.35           N  
ATOM    109  CA  SER A  65      32.130   8.131   4.315  1.00  8.84           C  
ATOM    110  C   SER A  65      30.817   7.862   3.590  1.00  8.60           C  
ATOM    111  O   SER A  65      29.878   8.664   3.639  1.00  8.27           O  
ATOM    112  CB  SER A  65      32.191   7.250   5.568  1.00  9.54           C  
ATOM    113  OG  SER A  65      31.018   7.423   6.342  1.00 10.43           O  
ATOM    114  N   GLU A  66      30.748   6.713   2.924  1.00  9.13           N  
ATOM    115  CA  GLU A  66      29.508   6.316   2.271  1.00 10.14           C  
ATOM    116  C   GLU A  66      28.356   6.255   3.260  1.00  9.39           C  
ATOM    117  O   GLU A  66      27.277   6.760   2.981  1.00  9.78           O  
ATOM    118  CB  GLU A  66      29.642   4.957   1.552  1.00 11.33           C  
ATOM    119  CG  GLU A  66      28.306   4.414   1.019  1.00 15.93           C  
ATOM    120  CD  GLU A  66      28.433   3.123   0.242  1.00 20.49           C  
ATOM    121  OE1 GLU A  66      28.933   2.141   0.821  1.00 23.61           O  
ATOM    122  OE2 GLU A  66      28.019   3.089  -0.939  1.00 25.25           O  
ATOM    123  N   SER A  67      28.585   5.632   4.407  1.00  9.14           N  
ATOM    124  CA  SER A  67      27.500   5.447   5.369  1.00 10.01           C  
ATOM    125  C   SER A  67      27.067   6.771   5.976  1.00  8.81           C  
ATOM    126  O   SER A  67      25.882   6.991   6.211  1.00  9.60           O  
ATOM    127  CB  SER A  67      27.874   4.458   6.473  1.00 11.61           C  
ATOM    128  OG  SER A  67      28.977   4.872   7.232  1.00 14.26           O  
ATOM    129  N   ARG A  68      28.017   7.667   6.221  1.00  8.58           N  
ATOM    130  CA  ARG A  68      27.693   8.974   6.785  1.00  8.32           C  
ATOM    131  C   ARG A  68      26.824   9.757   5.813  1.00  8.01           C  
ATOM    132  O   ARG A  68      25.796  10.325   6.203  1.00  7.98           O  
ATOM    133  CB  ARG A  68      28.955   9.785   7.107  1.00  8.65           C  
ATOM    134  CG  ARG A  68      28.649  11.128   7.782  1.00  8.97           C  
ATOM    135  CD  ARG A  68      28.366  10.995   9.266  1.00  9.86           C  
ATOM    136  NE  ARG A  68      29.617  10.789   9.969  1.00  9.23           N  
ATOM    137  CZ  ARG A  68      30.474  11.757  10.265  1.00  8.88           C  
ATOM    138  NH1 ARG A  68      30.148  13.034  10.076  1.00  9.14           N  
ATOM    139  NH2 ARG A  68      31.642  11.441  10.803  1.00  9.93           N  
ATOM    140  N   ILE A  69      27.229   9.824   4.555  1.00  7.72           N  
ATOM    141  CA  ILE A  69      26.494  10.586   3.566  1.00  7.90           C  
ATOM    142  C   ILE A  69      25.093  10.004   3.380  1.00  8.35           C  
ATOM    143  O   ILE A  69      24.116  10.740   3.335  1.00  8.77           O  
ATOM    144  CB  ILE A  69      27.280  10.669   2.238  1.00  8.29           C  
ATOM    145  CG1 ILE A  69      28.568  11.485   2.415  1.00  8.01           C  
ATOM    146  CG2 ILE A  69      26.414  11.244   1.119  1.00  9.03           C  
ATOM    147  CD1 ILE A  69      28.371  12.938   2.863  1.00  8.71           C  
ATOM    148  N   GLN A  70      24.978   8.687   3.301  1.00  9.03           N  
ATOM    149  CA  GLN A  70      23.663   8.065   3.154  1.00  9.97           C  
ATOM    150  C   GLN A  70      22.760   8.408   4.337  1.00  9.86           C  
ATOM    151  O   GLN A  70      21.590   8.713   4.159  1.00 10.39           O  
ATOM    152  CB  GLN A  70      23.802   6.548   3.000  1.00 11.39           C  
ATOM    153  CG  GLN A  70      24.345   6.122   1.634  1.00 15.15           C  
ATOM    154  CD  GLN A  70      24.586   4.613   1.506  1.00 19.50           C  
ATOM    155  OE1 GLN A  70      24.649   4.085   0.390  1.00 25.55           O  
ATOM    156  NE2 GLN A  70      24.751   3.927   2.639  1.00 21.99           N  
ATOM    157  N   SER A  71      23.313   8.370   5.539  1.00  9.16           N  
ATOM    158  CA  SER A  71      22.538   8.672   6.737  1.00  9.41           C  
ATOM    159  C   SER A  71      22.139  10.142   6.811  1.00  8.82           C  
ATOM    160  O   SER A  71      21.038  10.469   7.249  1.00  9.42           O  
ATOM    161  CB  SER A  71      23.330   8.303   7.980  1.00 10.04           C  
ATOM    162  OG  SER A  71      23.470   6.912   8.117  1.00 14.22           O  
ATOM    163  N   LEU A  72      23.025  11.035   6.393  1.00  7.98           N  
ATOM    164  CA  LEU A  72      22.700  12.454   6.385  1.00  7.91           C  
ATOM    165  C   LEU A  72      21.559  12.750   5.425  1.00  8.47           C  
ATOM    166  O   LEU A  72      20.677  13.537   5.738  1.00  9.08           O  
ATOM    167  CB  LEU A  72      23.937  13.291   6.066  1.00  7.31           C  
ATOM    168  CG  LEU A  72      24.961  13.356   7.203  1.00  7.27           C  
ATOM    169  CD1 LEU A  72      26.247  13.996   6.730  1.00  7.32           C  
ATOM    170  CD2 LEU A  72      24.421  14.102   8.419  1.00  7.53           C  
ATOM    171  N   HIS A  73      21.531  12.080   4.287  1.00  9.11           N  
ATOM    172  CA  HIS A  73      20.392  12.239   3.387  1.00 10.66           C  
ATOM    173  C   HIS A  73      19.112  11.631   3.938  1.00 11.06           C  
ATOM    174  O   HIS A  73      18.081  12.290   3.955  1.00 11.72           O  
ATOM    175  CB  HIS A  73      20.711  11.686   2.016  1.00 11.39           C  
ATOM    176  CG  HIS A  73      21.579  12.598   1.224  1.00 12.21           C  
ATOM    177  ND1 HIS A  73      21.077  13.656   0.500  1.00 15.54           N  
ATOM    178  CD2 HIS A  73      22.925  12.688   1.132  1.00 15.34           C  
ATOM    179  CE1 HIS A  73      22.073  14.308  -0.072  1.00 15.27           C  
ATOM    180  NE2 HIS A  73      23.205  13.744   0.302  1.00 16.00           N  
ATOM    181  N   LYS A  74      19.190  10.402   4.427  1.00 11.31           N  
ATOM    182  CA  LYS A  74      17.993   9.689   4.878  1.00 12.24           C  
ATOM    183  C   LYS A  74      17.318  10.385   6.052  1.00 11.98           C  
ATOM    184  O   LYS A  74      16.084  10.418   6.144  1.00 13.11           O  
ATOM    185  CB  LYS A  74      18.335   8.242   5.246  1.00 12.91           C  
ATOM    186  CG  LYS A  74      17.127   7.395   5.620  1.00 17.37           C  
ATOM    187  CD  LYS A  74      16.113   7.287   4.476  1.00 22.78           C  
ATOM    188  CE  LYS A  74      15.121   6.149   4.706  1.00 25.42           C  
ATOM    189  NZ  LYS A  74      14.309   5.856   3.497  1.00 28.28           N  
ATOM    190  N   SER A  75      18.123  10.938   6.951  1.00 11.07           N  
ATOM    191  CA  SER A  75      17.621  11.627   8.128  1.00 11.13           C  
ATOM    192  C   SER A  75      16.975  12.967   7.807  1.00 10.90           C  
ATOM    193  O   SER A  75      16.318  13.538   8.663  1.00 12.48           O  
ATOM    194  CB  SER A  75      18.750  11.878   9.116  1.00 11.16           C  
ATOM    195  OG  SER A  75      19.719  12.744   8.551  1.00 10.33           O  
ATOM    196  N   GLY A  76      17.193  13.497   6.610  1.00 10.13           N  
ATOM    197  CA  GLY A  76      16.750  14.833   6.275  1.00 10.10           C  
ATOM    198  C   GLY A  76      17.766  15.916   6.584  1.00  9.05           C  
ATOM    199  O   GLY A  76      17.565  17.064   6.193  1.00 10.10           O  
ATOM    200  N   LEU A  77      18.864  15.570   7.261  1.00  8.52           N  
ATOM    201  CA  LEU A  77      19.818  16.590   7.668  1.00  8.03           C  
ATOM    202  C   LEU A  77      20.554  17.221   6.488  1.00  7.21           C  
ATOM    203  O   LEU A  77      20.877  18.401   6.541  1.00  7.89           O  
ATOM    204  CB  LEU A  77      20.815  16.020   8.660  1.00  8.06           C  
ATOM    205  CG  LEU A  77      20.235  15.696  10.031  1.00  9.09           C  
ATOM    206  CD1 LEU A  77      21.234  14.940  10.837  1.00 10.45           C  
ATOM    207  CD2 LEU A  77      19.791  16.934  10.762  1.00 11.72           C  
ATOM    208  N   ALA A  78      20.808  16.462   5.424  1.00  7.83           N  
ATOM    209  CA  ALA A  78      21.581  16.968   4.295  1.00  8.47           C  
ATOM    210  C   ALA A  78      20.884  18.141   3.614  1.00  8.53           C  
ATOM    211  O   ALA A  78      21.533  19.032   3.088  1.00  9.36           O  
ATOM    212  CB  ALA A  78      21.839  15.872   3.289  1.00  9.03           C  
ATOM    213  N   ALA A  79      19.550  18.145   3.630  1.00  9.47           N  
ATOM    214  CA  ALA A  79      18.752  19.184   2.970  1.00 10.09           C  
ATOM    215  C   ALA A  79      18.653  20.471   3.783  1.00  9.92           C  
ATOM    216  O   ALA A  79      18.143  21.473   3.292  1.00 11.33           O  
ATOM    217  CB  ALA A  79      17.360  18.646   2.708  1.00 10.82           C  
ATOM    218  N   GLN A  80      19.116  20.441   5.033  1.00  9.26           N  
ATOM    219  CA  GLN A  80      19.085  21.596   5.919  1.00  9.29           C  
ATOM    220  C   GLN A  80      20.326  22.465   5.664  1.00  9.11           C  
ATOM    221  O   GLN A  80      20.956  22.333   4.625  1.00  9.99           O  
ATOM    222  CB  GLN A  80      18.933  21.118   7.359  1.00  9.29           C  
ATOM    223  CG  GLN A  80      17.629  20.324   7.551  1.00 10.26           C  
ATOM    224  CD  GLN A  80      17.451  19.716   8.934  1.00 12.16           C  
ATOM    225  OE1 GLN A  80      18.187  20.005   9.855  1.00 12.95           O  
ATOM    226  NE2 GLN A  80      16.436  18.870   9.074  1.00 15.56           N  
ATOM    227  N   GLU A  81      20.643  23.377   6.573  1.00  8.73           N  
ATOM    228  CA  GLU A  81      21.690  24.368   6.333  1.00  8.29           C  
ATOM    229  C   GLU A  81      23.011  23.942   6.965  1.00  6.84           C  
ATOM    230  O   GLU A  81      23.063  23.575   8.143  1.00  6.87           O  
ATOM    231  CB  GLU A  81      21.292  25.723   6.915  1.00  9.67           C  
ATOM    232  CG  GLU A  81      20.281  26.580   6.168  1.00 13.76           C  
ATOM    233  CD  GLU A  81      20.492  28.042   6.545  1.00 16.85           C  
ATOM    234  OE1 GLU A  81      21.004  28.826   5.719  1.00 18.18           O  
ATOM    235  OE2 GLU A  81      20.210  28.388   7.713  1.00 19.24           O  
ATOM    236  N   TRP A  82      24.055  24.069   6.162  1.00  6.78           N  
ATOM    237  CA  TRP A  82      25.428  23.727   6.498  1.00  5.64           C  
ATOM    238  C   TRP A  82      26.323  24.892   6.106  1.00  5.61           C  
ATOM    239  O   TRP A  82      25.987  25.707   5.245  1.00  6.42           O  
ATOM    240  CB  TRP A  82      25.861  22.476   5.715  1.00  6.41           C  
ATOM    241  CG  TRP A  82      25.084  21.243   6.054  1.00  6.10           C  
ATOM    242  CD1 TRP A  82      23.780  20.996   5.757  1.00  6.25           C  
ATOM    243  CD2 TRP A  82      25.558  20.090   6.748  1.00  5.81           C  
ATOM    244  NE1 TRP A  82      23.407  19.771   6.247  1.00  6.29           N  
ATOM    245  CE2 TRP A  82      24.477  19.195   6.865  1.00  6.66           C  
ATOM    246  CE3 TRP A  82      26.794  19.722   7.294  1.00  6.08           C  
ATOM    247  CZ2 TRP A  82      24.599  17.961   7.483  1.00  7.04           C  
ATOM    248  CZ3 TRP A  82      26.914  18.497   7.903  1.00  6.38           C  
ATOM    249  CH2 TRP A  82      25.832  17.632   7.997  1.00  7.28           C  
ATOM    250  N   VAL A  83      27.495  24.938   6.725  1.00  5.50           N  
ATOM    251  CA  VAL A  83      28.558  25.853   6.334  1.00  5.35           C  
ATOM    252  C   VAL A  83      29.852  25.080   6.127  1.00  5.46           C  
ATOM    253  O   VAL A  83      30.092  24.049   6.745  1.00  6.18           O  
ATOM    254  CB  VAL A  83      28.792  26.978   7.383  1.00  5.76           C  
ATOM    255  CG1 VAL A  83      27.590  27.932   7.390  1.00  6.61           C  
ATOM    256  CG2 VAL A  83      29.071  26.406   8.787  1.00  6.54           C  
ATOM    257  N   ALA A  84      30.707  25.624   5.272  1.00  5.43           N  
ATOM    258  CA  ALA A  84      32.110  25.240   5.201  1.00  5.06           C  
ATOM    259  C   ALA A  84      32.946  26.392   5.715  1.00  5.11           C  
ATOM    260  O   ALA A  84      32.764  27.514   5.253  1.00  5.97           O  
ATOM    261  CB  ALA A  84      32.527  24.920   3.791  1.00  6.31           C  
ATOM    262  N   CYS A  85      33.852  26.124   6.654  1.00  4.78           N  
ATOM    263  CA  CYS A  85      34.795  27.114   7.154  1.00  4.73           C  
ATOM    264  C   CYS A  85      36.205  26.576   6.987  1.00  4.63           C  
ATOM    265  O   CYS A  85      36.421  25.390   6.800  1.00  4.76           O  
ATOM    266  CB  CYS A  85      34.535  27.427   8.627  1.00  5.52           C  
ATOM    267  SG  CYS A  85      32.869  28.017   9.012  1.00  7.11           S  
ATOM    268  N   GLU A  86      37.190  27.455   7.078  1.00  4.62           N  
ATOM    269  CA  GLU A  86      38.571  27.044   6.878  1.00  4.45           C  
ATOM    270  C   GLU A  86      39.020  26.076   7.969  1.00  4.37           C  
ATOM    271  O   GLU A  86      38.754  26.297   9.154  1.00  4.85           O  
ATOM    272  CB  GLU A  86      39.465  28.281   6.845  1.00  5.19           C  
ATOM    273  CG  GLU A  86      40.931  27.959   6.607  1.00  4.98           C  
ATOM    274  CD  GLU A  86      41.812  29.159   6.344  1.00  6.64           C  
ATOM    275  OE1 GLU A  86      41.342  30.314   6.471  1.00  7.94           O  
ATOM    276  OE2 GLU A  86      43.012  28.933   6.041  1.00  6.62           O  
ATOM    277  N   LYS A  87      39.723  25.026   7.552  1.00  3.83           N  
ATOM    278  CA  LYS A  87      40.425  24.140   8.472  1.00  3.94           C  
ATOM    279  C   LYS A  87      41.850  24.659   8.607  1.00  3.95           C  
ATOM    280  O   LYS A  87      42.583  24.714   7.619  1.00  4.64           O  
ATOM    281  CB  LYS A  87      40.406  22.700   7.983  1.00  4.21           C  
ATOM    282  CG  LYS A  87      40.812  21.716   9.060  1.00  4.82           C  
ATOM    283  CD  LYS A  87      40.599  20.285   8.636  1.00  5.64           C  
ATOM    284  CE  LYS A  87      40.745  19.292   9.788  1.00  6.04           C  
ATOM    285  NZ  LYS A  87      42.137  19.250  10.300  1.00  6.11           N  
ATOM    286  N   VAL A  88      42.191  25.106   9.813  1.00  4.11           N  
ATOM    287  CA  VAL A  88      43.484  25.722  10.086  1.00  4.45           C  
ATOM    288  C   VAL A  88      44.463  24.675  10.618  1.00  3.93           C  
ATOM    289  O   VAL A  88      44.131  23.888  11.506  1.00  4.74           O  
ATOM    290  CB  VAL A  88      43.302  26.914  11.068  1.00  4.47           C  
ATOM    291  CG1 VAL A  88      44.644  27.493  11.514  1.00  5.48           C  
ATOM    292  CG2 VAL A  88      42.459  27.975  10.373  1.00  5.25           C  
ATOM    293  N   HIS A  89      45.670  24.687  10.054  1.00  4.36           N  
ATOM    294  CA  HIS A  89      46.730  23.760  10.416  1.00  4.44           C  
ATOM    295  C   HIS A  89      47.626  24.352  11.490  1.00  4.47           C  
ATOM    296  O   HIS A  89      48.645  24.992  11.200  1.00  5.45           O  
ATOM    297  CB  HIS A  89      47.539  23.442   9.194  1.00  5.39           C  
ATOM    298  CG  HIS A  89      48.572  22.385   9.383  1.00  5.84           C  
ATOM    299  ND1 HIS A  89      48.273  21.121   9.793  1.00  7.90           N  
ATOM    300  CD2 HIS A  89      49.912  22.427   9.275  1.00  9.92           C  
ATOM    301  CE1 HIS A  89      49.367  20.386   9.781  1.00  7.28           C  
ATOM    302  NE2 HIS A  89      50.363  21.139   9.389  1.00  7.99           N  
ATOM    303  N   GLY A  90      47.212  24.171  12.734  1.00  4.84           N  
ATOM    304  CA  GLY A  90      47.977  24.569  13.905  1.00  4.82           C  
ATOM    305  C   GLY A  90      47.968  23.416  14.884  1.00  4.31           C  
ATOM    306  O   GLY A  90      48.231  22.268  14.522  1.00  5.30           O  
ATOM    307  N   THR A  91      47.650  23.727  16.133  1.00  4.47           N  
ATOM    308  CA  THR A  91      47.488  22.698  17.144  1.00  5.08           C  
ATOM    309  C   THR A  91      46.227  23.006  17.947  1.00  4.89           C  
ATOM    310  O   THR A  91      45.848  24.161  18.123  1.00  4.88           O  
ATOM    311  CB  THR A  91      48.776  22.587  17.981  1.00  5.80           C  
ATOM    312  OG1 THR A  91      48.738  21.401  18.776  1.00  7.22           O  
ATOM    313  CG2 THR A  91      48.947  23.756  18.929  1.00  7.35           C  
ATOM    314  N   ASN A  92      45.560  21.960  18.391  1.00  4.91           N  
ATOM    315  CA  ASN A  92      44.301  22.113  19.093  1.00  5.14           C  
ATOM    316  C   ASN A  92      44.489  22.853  20.414  1.00  4.84           C  
ATOM    317  O   ASN A  92      45.438  22.578  21.157  1.00  5.79           O  
ATOM    318  CB  ASN A  92      43.702  20.727  19.353  1.00  5.44           C  
ATOM    319  CG  ASN A  92      42.298  20.794  19.874  1.00  6.04           C  
ATOM    320  OD1 ASN A  92      42.067  20.782  21.077  1.00  9.32           O  
ATOM    321  ND2 ASN A  92      41.350  20.907  18.972  1.00  6.71           N  
ATOM    322  N   PHE A  93      43.585  23.765  20.710  1.00  5.58           N  
ATOM    323  CA  PHE A  93      43.744  24.595  21.888  1.00  6.28           C  
ATOM    324  C   PHE A  93      42.372  24.872  22.491  1.00  6.64           C  
ATOM    325  O   PHE A  93      41.336  24.938  21.797  1.00 11.84           O  
ATOM    326  CB  PHE A  93      44.514  25.854  21.484  1.00  6.74           C  
ATOM    327  CG  PHE A  93      45.044  26.656  22.642  1.00  6.65           C  
ATOM    328  CD1 PHE A  93      46.235  26.319  23.260  1.00  8.52           C  
ATOM    329  CD2 PHE A  93      44.360  27.767  23.113  1.00  6.80           C  
ATOM    330  CE1 PHE A  93      46.720  27.067  24.325  1.00  9.10           C  
ATOM    331  CE2 PHE A  93      44.847  28.519  24.178  1.00  8.18           C  
ATOM    332  CZ  PHE A  93      46.031  28.162  24.776  1.00  8.89           C  
ATOM    333  N   GLY A  94      42.300  24.959  23.791  1.00  5.27           N  
ATOM    334  CA  GLY A  94      41.063  25.306  24.447  1.00  5.53           C  
ATOM    335  C   GLY A  94      41.307  26.336  25.528  1.00  5.32           C  
ATOM    336  O   GLY A  94      42.275  26.237  26.289  1.00  6.37           O  
ATOM    337  N   ILE A  95      40.417  27.313  25.600  1.00  5.38           N  
ATOM    338  CA  ILE A  95      40.413  28.316  26.652  1.00  5.67           C  
ATOM    339  C   ILE A  95      39.185  28.061  27.513  1.00  6.00           C  
ATOM    340  O   ILE A  95      38.057  28.034  27.030  1.00  6.51           O  
ATOM    341  CB  ILE A  95      40.343  29.725  26.042  1.00  6.28           C  
ATOM    342  CG1 ILE A  95      41.534  29.973  25.112  1.00  7.15           C  
ATOM    343  CG2 ILE A  95      40.270  30.771  27.132  1.00  7.90           C  
ATOM    344  CD1 ILE A  95      41.367  31.118  24.145  1.00  8.23           C  
ATOM    345  N   TYR A  96      39.433  27.859  28.798  1.00  6.44           N  
ATOM    346  CA  TYR A  96      38.417  27.482  29.757  1.00  7.18           C  
ATOM    347  C   TYR A  96      38.202  28.577  30.791  1.00  7.54           C  
ATOM    348  O   TYR A  96      39.158  29.130  31.312  1.00  8.77           O  
ATOM    349  CB  TYR A  96      38.872  26.235  30.517  1.00  7.82           C  
ATOM    350  CG  TYR A  96      38.843  24.964  29.703  1.00  7.38           C  
ATOM    351  CD1 TYR A  96      37.934  23.971  29.993  1.00  7.62           C  
ATOM    352  CD2 TYR A  96      39.720  24.749  28.639  1.00  7.58           C  
ATOM    353  CE1 TYR A  96      37.871  22.814  29.270  1.00  8.09           C  
ATOM    354  CE2 TYR A  96      39.659  23.587  27.890  1.00  7.82           C  
ATOM    355  CZ  TYR A  96      38.726  22.620  28.218  1.00  7.05           C  
ATOM    356  OH  TYR A  96      38.628  21.447  27.518  1.00  8.76           O  
ATOM    357  N   LEU A  97      36.940  28.832  31.123  1.00  7.39           N  
ATOM    358  CA  LEU A  97      36.602  29.617  32.301  1.00  7.48           C  
ATOM    359  C   LEU A  97      35.731  28.717  33.170  1.00  8.25           C  
ATOM    360  O   LEU A  97      34.679  28.249  32.754  1.00  9.30           O  
ATOM    361  CB  LEU A  97      35.855  30.893  31.933  1.00  8.20           C  
ATOM    362  CG  LEU A  97      35.619  31.870  33.083  1.00 10.44           C  
ATOM    363  CD1 LEU A  97      36.893  32.376  33.661  1.00 12.92           C  
ATOM    364  CD2 LEU A  97      34.813  33.035  32.587  1.00 13.72           C  
ATOM    365  N   ILE A  98      36.192  28.453  34.388  1.00  9.59           N  
ATOM    366  CA  ILE A  98      35.572  27.483  35.291  1.00 11.31           C  
ATOM    367  C   ILE A  98      35.063  28.195  36.537  1.00 13.03           C  
ATOM    368  O   ILE A  98      35.846  28.803  37.223  1.00 13.62           O  
ATOM    369  CB  ILE A  98      36.626  26.432  35.705  1.00 11.51           C  
ATOM    370  CG1 ILE A  98      37.266  25.760  34.478  1.00 14.23           C  
ATOM    371  CG2 ILE A  98      36.041  25.420  36.713  1.00 12.82           C  
ATOM    372  CD1 ILE A  98      36.407  24.977  33.692  1.00 17.84           C  
ATOM    373  N   ASN A  99      33.755  28.126  36.771  1.00 13.87           N  
ATOM    374  CA  ASN A  99      33.124  28.671  37.963  1.00 15.86           C  
ATOM    375  C   ASN A  99      33.287  27.740  39.147  1.00 16.91           C  
ATOM    376  O   ASN A  99      33.064  26.526  39.060  1.00 18.08           O  
ATOM    377  CB  ASN A  99      31.634  28.888  37.739  1.00 16.88           C  
ATOM    378  CG  ASN A  99      30.971  29.502  38.937  1.00 18.71           C  
ATOM    379  OD1 ASN A  99      30.409  28.800  39.777  1.00 21.30           O  
ATOM    380  ND2 ASN A  99      31.073  30.813  39.052  1.00 20.44           N  
ATOM    381  N   GLN A 100      33.694  28.340  40.249  1.00 17.83           N  
ATOM    382  CA  GLN A 100      33.844  27.655  41.518  1.00 19.07           C  
ATOM    383  C   GLN A 100      33.106  28.467  42.579  1.00 19.98           C  
ATOM    384  O   GLN A 100      33.696  28.935  43.547  1.00 20.42           O  
ATOM    385  CB  GLN A 100      35.326  27.532  41.833  1.00 19.00           C  
ATOM    386  CG  GLN A 100      36.072  26.694  40.803  1.00 20.95           C  
ATOM    387  CD  GLN A 100      37.569  26.712  40.988  1.00 23.06           C  
ATOM    388  OE1 GLN A 100      38.112  27.582  41.673  1.00 25.81           O  
ATOM    389  NE2 GLN A 100      38.248  25.753  40.365  1.00 24.14           N  
ATOM    390  N   GLY A 101      31.803  28.633  42.384  1.00 21.10           N  
ATOM    391  CA  GLY A 101      30.973  29.381  43.319  1.00 21.76           C  
ATOM    392  C   GLY A 101      31.099  30.883  43.122  1.00 22.45           C  
ATOM    393  O   GLY A 101      30.743  31.399  42.061  1.00 22.91           O  
ATOM    394  N   ASP A 102      31.617  31.594  44.126  1.00 23.19           N  
ATOM    395  CA  ASP A 102      31.788  33.054  44.029  1.00 23.87           C  
ATOM    396  C   ASP A 102      33.144  33.431  43.435  1.00 23.12           C  
ATOM    397  O   ASP A 102      33.490  34.611  43.376  1.00 23.97           O  
ATOM    398  CB  ASP A 102      31.616  33.777  45.382  1.00 24.59           C  
ATOM    399  CG  ASP A 102      31.047  32.897  46.460  1.00 27.28           C  
ATOM    400  OD1 ASP A 102      29.909  32.403  46.292  1.00 31.96           O  
ATOM    401  OD2 ASP A 102      31.670  32.656  47.515  1.00 30.79           O  
"""

    def pdbqt_file_str(self):
        """Get a PDBQT-formatted string.
        
            Returns:
                A string, PDBQT formatted.
        """

        return """ATOM      1  N   PRO     1     111.971  33.330  97.856  1.00  0.00    -0.062 N 
ATOM      2  HN1 PRO     1     111.233  33.160  97.188  1.00  0.00     0.278 HD
ATOM      3  HN2 PRO     1     111.849  33.063  98.822  1.00  0.00     0.278 HD
ATOM      4  CD  PRO     1     113.196  32.635  97.440  1.00  0.00     0.233 C 
ATOM      5  CG  PRO     1     114.270  33.722  97.624  1.00  0.00     0.030 C 
ATOM      6  CB  PRO     1     113.573  35.036  97.237  1.00  0.00     0.044 C 
ATOM      7  CA  PRO     1     112.160  34.798  97.760  1.00  0.00     0.277 C 
ATOM      8  C   PRO     1     111.127  35.407  96.816  1.00  0.00     0.249 C 
ATOM      9  O   PRO     1     110.966  34.956  95.677  1.00  0.00    -0.271 OA
ATOM     10  N   ARG     2     110.425  36.428  97.295  1.00  0.00    -0.346 N 
ATOM     11  HN  ARG     2     110.557  36.728  98.250  1.00  0.00     0.163 HD
ATOM     12  CA  ARG     2     109.409  37.101  96.493  1.00  0.00     0.176 C 
ATOM     13  CB  ARG     2     108.022  36.629  96.916  1.00  0.00     0.036 C 
ATOM     14  CG  ARG     2     107.636  35.304  96.297  1.00  0.00     0.023 C 
ATOM     15  CD  ARG     2     106.771  34.492  97.234  1.00  0.00     0.138 C 
ATOM     16  NE  ARG     2     105.707  35.291  97.831  1.00  0.00    -0.227 N 
ATOM     17  HE  ARG     2     105.648  36.263  97.564  1.00  0.00     0.177 HD
ATOM     18  CZ  ARG     2     104.814  34.815  98.691  1.00  0.00     0.665 C 
ATOM     19  NH1 ARG     2     104.857  33.539  99.055  1.00  0.00    -0.235 N 
ATOM     20 1HH1 ARG     2     105.570  32.933  98.676  1.00  0.00     0.174 HD
ATOM     21 2HH1 ARG     2     104.176  33.180  99.709  1.00  0.00     0.174 HD
ATOM     22  NH2 ARG     2     103.884  35.612  99.196  1.00  0.00    -0.235 N 
ATOM     23 1HH2 ARG     2     103.855  36.585  98.925  1.00  0.00     0.174 HD
ATOM     24 2HH2 ARG     2     103.208  35.243  99.850  1.00  0.00     0.174 HD
ATOM     25  C   ARG     2     109.516  38.616  96.608  1.00  0.00     0.241 C 
ATOM     26  O   ARG     2     108.956  39.229  97.515  1.00  0.00    -0.271 OA
ATOM     27  N   ILE     3     110.241  39.214  95.669  1.00  0.00    -0.346 N 
ATOM     28  HN  ILE     3     110.670  38.667  94.936  1.00  0.00     0.163 HD
ATOM     29  CA  ILE     3     110.434  40.657  95.656  1.00  0.00     0.180 C 
ATOM     30  CB  ILE     3     111.680  41.047  94.840  1.00  0.00     0.013 C 
ATOM     31  CG2 ILE     3     112.019  42.517  95.066  1.00  0.00     0.012 C 
ATOM     32  CG1 ILE     3     112.865  40.174  95.240  1.00  0.00     0.002 C 
ATOM     33  CD1 ILE     3     114.048  40.287  94.287  1.00  0.00     0.005 C 
ATOM     34  C   ILE     3     109.248  41.380  95.031  1.00  0.00     0.241 C 
ATOM     35  O   ILE     3     108.776  40.998  93.954  1.00  0.00    -0.271 OA
ATOM     36  N   LYS     4     108.778  42.425  95.710  1.00  0.00    -0.346 N 
ATOM     37  HN  LYS     4     109.206  42.676  96.590  1.00  0.00     0.163 HD
ATOM     38  CA  LYS     4     107.674  43.243  95.220  1.00  0.00     0.176 C 
ATOM     39  CB  LYS     4     106.541  43.297  96.246  1.00  0.00     0.035 C 
ATOM     40  CG  LYS     4     105.366  44.203  95.858  1.00  0.00     0.004 C 
ATOM     41  CD  LYS     4     104.603  43.689  94.631  1.00  0.00     0.027 C 
ATOM     42  CE  LYS     4     103.394  44.567  94.318  1.00  0.00     0.229 C 
ATOM     43  NZ  LYS     4     102.627  44.109  93.118  1.00  0.00    -0.079 N 
ATOM     44  HZ1 LYS     4     102.290  43.170  93.272  1.00  0.00     0.274 HD
ATOM     45  HZ2 LYS     4     101.843  44.727  92.963  1.00  0.00     0.274 HD
ATOM     46  HZ3 LYS     4     103.230  44.121  92.308  1.00  0.00     0.274 HD
ATOM     47  C   LYS     4     108.240  44.640  95.010  1.00  0.00     0.241 C 
ATOM     48  O   LYS     4     109.013  45.123  95.833  1.00  0.00    -0.271 OA
ATOM     49  N   LYS     5     107.877  45.282  93.906  1.00  0.00    -0.346 N 
ATOM     50  HN  LYS     5     107.282  44.842  93.219  1.00  0.00     0.163 HD
ATOM     51  CA  LYS     5     108.371  46.627  93.634  1.00  0.00     0.176 C 
ATOM     52  CB  LYS     5     108.612  46.828  92.131  1.00  0.00     0.035 C 
ATOM     53  CG  LYS     5     110.049  46.571  91.693  1.00  0.00     0.004 C 
ATOM     54  CD  LYS     5     110.323  47.140  90.308  1.00  0.00     0.027 C 
ATOM     55  CE  LYS     5     111.799  47.019  89.941  1.00  0.00     0.229 C 
ATOM     56  NZ  LYS     5     112.121  47.763  88.688  1.00  0.00    -0.079 N 
ATOM     57  HZ1 LYS     5     111.904  48.742  88.812  1.00  0.00     0.274 HD
ATOM     58  HZ2 LYS     5     113.104  47.660  88.478  1.00  0.00     0.274 HD
ATOM     59  HZ3 LYS     5     111.573  47.391  87.925  1.00  0.00     0.274 HD
ATOM     60  C   LYS     5     107.404  47.689  94.138  1.00  0.00     0.241 C 
ATOM     61  O   LYS     5     106.189  47.536  94.024  1.00  0.00    -0.271 OA
ATOM     62  N   PHE     6     107.949  48.762  94.699  1.00  0.00    -0.346 N 
ATOM     63  HN  PHE     6     108.950  48.818  94.825  1.00  0.00     0.163 HD
ATOM     64  CA  PHE     6     107.128  49.850  95.208  1.00  0.00     0.180 C 
ATOM     65  CB  PHE     6     107.113  49.864  96.745  1.00  0.00     0.073 C 
ATOM     66  CG  PHE     6     106.232  48.808  97.351  1.00  0.00    -0.056 A 
ATOM     67  CD1 PHE     6     106.777  47.630  97.854  1.00  0.00     0.007 A 
ATOM     68  CE1 PHE     6     105.955  46.617  98.370  1.00  0.00     0.001 A 
ATOM     69  CZ  PHE     6     104.568  46.787  98.385  1.00  0.00     0.000 A 
ATOM     70  CE2 PHE     6     104.013  47.969  97.885  1.00  0.00     0.001 A 
ATOM     71  CD2 PHE     6     104.848  48.971  97.372  1.00  0.00     0.007 A 
ATOM     72  C   PHE     6     107.563  51.212  94.706  1.00  0.00     0.241 C 
ATOM     73  O   PHE     6     108.664  51.684  95.000  1.00  0.00    -0.271 OA
ATOM     74  N   ALA     7     106.689  51.840  93.932  1.00  0.00    -0.346 N 
ATOM     75  HN  ALA     7     105.802  51.411  93.708  1.00  0.00     0.163 HD
ATOM     76  CA  ALA     7     106.960  53.172  93.418  1.00  0.00     0.172 C 
ATOM     77  CB  ALA     7     106.395  53.323  92.018  1.00  0.00     0.042 C 
ATOM     78  C   ALA     7     106.226  54.085  94.387  1.00  0.00     0.240 C 
ATOM     79  O   ALA     7     105.016  53.947  94.572  1.00  0.00    -0.271 OA
ATOM     80  N   ILE     8     106.956  55.005  95.011  1.00  0.00    -0.346 N 
ATOM     81  HN  ILE     8     107.951  55.068  94.849  1.00  0.00     0.163 HD
ATOM     82  CA  ILE     8     106.355  55.911  95.978  1.00  0.00     0.180 C 
ATOM     83  CB  ILE     8     106.786  55.552  97.408  1.00  0.00     0.013 C 
ATOM     84  CG2 ILE     8     106.246  56.593  98.390  1.00  0.00     0.012 C 
ATOM     85  CG1 ILE     8     106.293  54.144  97.761  1.00  0.00     0.002 C 
ATOM     86  CD1 ILE     8     106.691  53.663  99.144  1.00  0.00     0.005 C 
ATOM     87  C   ILE     8     106.678  57.376  95.749  1.00  0.00     0.241 C 
ATOM     88  O   ILE     8     107.829  57.745  95.534  1.00  0.00    -0.271 OA
ATOM     89  N   TYR     9     105.648  58.206  95.813  1.00  0.00    -0.346 N 
ATOM     90  HN  TYR     9     104.727  57.826  95.977  1.00  0.00     0.163 HD
ATOM     91  CA  TYR     9     105.779  59.640  95.641  1.00  0.00     0.180 C 
ATOM     92  CB  TYR     9     104.377  60.256  95.700  1.00  0.00     0.073 C 
ATOM     93  CG  TYR     9     104.312  61.765  95.625  1.00  0.00    -0.056 A 
ATOM     94  CD1 TYR     9     104.562  62.555  96.745  1.00  0.00     0.010 A 
ATOM     95  CE1 TYR     9     104.460  63.945  96.684  1.00  0.00     0.037 A 
ATOM     96  CZ  TYR     9     104.108  64.557  95.490  1.00  0.00     0.065 A 
ATOM     97  OH  TYR     9     103.993  65.930  95.406  1.00  0.00    -0.361 OA
ATOM     98  HH  TYR     9     104.179  66.374  96.237  1.00  0.00     0.217 HD
ATOM     99  CE2 TYR     9     103.861  63.792  94.366  1.00  0.00     0.037 A 
ATOM    100  CD2 TYR     9     103.964  62.403  94.438  1.00  0.00     0.010 A 
ATOM    101  C   TYR     9     106.672  60.201  96.754  1.00  0.00     0.241 C 
ATOM    102  O   TYR     9     106.666  59.681  97.876  1.00  0.00    -0.271 OA
ATOM    103  N   ARG    10     107.436  61.253  96.451  1.00  0.00    -0.346 N 
ATOM    104  HN  ARG    10     107.431  61.632  95.515  1.00  0.00     0.163 HD
ATOM    105  CA  ARG    10     108.323  61.869  97.446  1.00  0.00     0.176 C 
ATOM    106  CB  ARG    10     109.714  61.215  97.411  1.00  0.00     0.036 C 
ATOM    107  CG  ARG    10     109.741  59.747  97.773  1.00  0.00     0.023 C 
ATOM    108  CD  ARG    10     109.302  59.513  99.207  1.00  0.00     0.138 C 
ATOM    109  NE  ARG    10     110.419  59.203 100.095  1.00  0.00    -0.227 N 
ATOM    110  HE  ARG    10     110.558  58.230 100.328  1.00  0.00     0.177 HD
ATOM    111  CZ  ARG    10     111.259  60.101 100.597  1.00  0.00     0.665 C 
ATOM    112  NH1 ARG    10     111.126  61.390 100.309  1.00  0.00    -0.235 N 
ATOM    113 1HH1 ARG    10     110.379  61.694  99.702  1.00  0.00     0.174 HD
ATOM    114 2HH1 ARG    10     111.773  62.061 100.699  1.00  0.00     0.174 HD
ATOM    115  NH2 ARG    10     112.229  59.708 101.403  1.00  0.00    -0.235 N 
ATOM    116 1HH2 ARG    10     112.326  58.729 101.632  1.00  0.00     0.174 HD
ATOM    117 2HH2 ARG    10     112.870  60.388 101.786  1.00  0.00     0.174 HD
ATOM    118  C   ARG    10     108.500  63.369  97.242  1.00  0.00     0.241 C 
ATOM    119  O   ARG    10     108.377  63.866  96.134  1.00  0.00    -0.271 OA
ATOM    120  N   TRP    11     108.779  64.086  98.325  1.00  0.00    -0.346 N 
ATOM    121  HN  TRP    11     108.798  63.641  99.231  1.00  0.00     0.163 HD
ATOM    122  CA  TRP    11     109.021  65.527  98.255  1.00  0.00     0.181 C 
ATOM    123  CB  TRP    11     107.807  66.308  97.786  1.00  0.00     0.075 C 
ATOM    124  CG  TRP    11     108.163  67.754  97.699  1.00  0.00    -0.028 A 
ATOM    125  CD1 TRP    11     109.052  68.322  96.832  1.00  0.00     0.096 A 
ATOM    126  NE1 TRP    11     109.180  69.664  97.089  1.00  0.00    -0.365 N 
ATOM    127  HE1 TRP    11     109.776  70.304  96.584  1.00  0.00     0.165 HD
ATOM    128  CE2 TRP    11     108.368  69.995  98.139  1.00  0.00     0.042 A 
ATOM    129  CZ2 TRP    11     108.156  71.232  98.762  1.00  0.00     0.030 A 
ATOM    130  CH2 TRP    11     107.269  71.274  99.810  1.00  0.00     0.002 A 
ATOM    131  CZ3 TRP    11     106.593  70.116 100.249  1.00  0.00     0.001 A 
ATOM    132  CE3 TRP    11     106.805  68.879  99.626  1.00  0.00     0.014 A 
ATOM    133  CD2 TRP    11     107.706  68.812  98.552  1.00  0.00    -0.002 A 
ATOM    134  C   TRP    11     109.453  66.109  99.588  1.00  0.00     0.241 C 
ATOM    135  O   TRP    11     108.675  66.165 100.542  1.00  0.00    -0.271 OA
ATOM    136  N   ASP    12     110.692  66.578  99.633  1.00  0.00    -0.345 N 
ATOM    137  HN  ASP    12     111.290  66.503  98.823  1.00  0.00     0.163 HD
ATOM    138  CA  ASP    12     111.246  67.136 100.844  1.00  0.00     0.187 C 
ATOM    139  CB  ASP    12     112.696  66.660 100.996  1.00  0.00     0.147 C 
ATOM    140  CG  ASP    12     113.273  66.962 102.366  1.00  0.00     0.175 C 
ATOM    141  OD1 ASP    12     114.264  66.303 102.753  1.00  0.00    -0.648 OA
ATOM    142  OD2 ASP    12     112.745  67.863 103.053  1.00  0.00    -0.648 OA
ATOM    143  C   ASP    12     111.167  68.660 100.870  1.00  0.00     0.243 C 
ATOM    144  O   ASP    12     111.764  69.341 100.044  1.00  0.00    -0.271 OA
ATOM    145  N   PRO    13     110.411  69.212 101.829  1.00  0.00    -0.337 N 
ATOM    146  CD  PRO    13     109.584  68.479 102.803  1.00  0.00     0.127 C 
ATOM    147  CG  PRO    13     108.483  69.456 103.063  1.00  0.00     0.022 C 
ATOM    148  CB  PRO    13     109.192  70.780 103.088  1.00  0.00     0.037 C 
ATOM    149  CA  PRO    13     110.244  70.663 101.982  1.00  0.00     0.179 C 
ATOM    150  C   PRO    13     111.560  71.317 102.403  1.00  0.00     0.241 C 
ATOM    151  O   PRO    13     111.801  72.492 102.123  1.00  0.00    -0.271 OA
ATOM    152  N   ASP    14     112.402  70.544 103.087  1.00  0.00    -0.346 N 
ATOM    153  HN  ASP    14     112.142  69.594 103.310  1.00  0.00     0.163 HD
ATOM    154  CA  ASP    14     113.685  71.043 103.575  1.00  0.00     0.186 C 
ATOM    155  CB  ASP    14     113.973  70.530 104.993  1.00  0.00     0.147 C 
ATOM    156  CG  ASP    14     112.914  70.937 105.995  1.00  0.00     0.175 C 
ATOM    157  OD1 ASP    14     111.930  70.184 106.158  1.00  0.00    -0.648 OA
ATOM    158  OD2 ASP    14     113.062  72.012 106.614  1.00  0.00    -0.648 OA
ATOM    159  C   ASP    14     114.857  70.668 102.684  1.00  0.00     0.241 C 
ATOM    160  O   ASP    14     115.967  70.453 103.173  1.00  0.00    -0.271 OA
ATOM    161  N   LYS    15     114.610  70.571 101.383  1.00  0.00    -0.346 N 
ATOM    162  HN  LYS    15     113.673  70.705 101.032  1.00  0.00     0.163 HD
ATOM    163  CA  LYS    15     115.665  70.246 100.430  1.00  0.00     0.176 C 
ATOM    164  CB  LYS    15     115.567  68.787  99.962  1.00  0.00     0.035 C 
ATOM    165  CG  LYS    15     116.768  68.343  99.140  1.00  0.00     0.004 C 
ATOM    166  CD  LYS    15     116.692  66.881  98.730  1.00  0.00     0.027 C 
ATOM    167  CE  LYS    15     117.875  66.516  97.834  1.00  0.00     0.229 C 
ATOM    168  NZ  LYS    15     117.775  65.146  97.249  1.00  0.00    -0.079 N 
ATOM    169  HZ1 LYS    15     116.938  65.080  96.688  1.00  0.00     0.274 HD
ATOM    170  HZ2 LYS    15     118.583  64.967  96.670  1.00  0.00     0.274 HD
ATOM    171  HZ3 LYS    15     117.736  64.464  97.993  1.00  0.00     0.274 HD
ATOM    172  C   LYS    15     115.528  71.191  99.246  1.00  0.00     0.241 C 
ATOM    173  O   LYS    15     114.843  70.895  98.270  1.00  0.00    -0.271 OA
ATOM    174  N   THR    16     116.181  72.341  99.360  1.00  0.00    -0.344 N 
ATOM    175  HN  THR    16     116.695  72.527 100.209  1.00  0.00     0.163 HD
ATOM    176  CA  THR    16     116.146  73.375  98.333  1.00  0.00     0.205 C 
ATOM    177  CB  THR    16     117.419  74.246  98.385  1.00  0.00     0.146 C 
ATOM    178  CG2 THR    16     117.226  75.511  97.556  1.00  0.00     0.042 C 
ATOM    179  OG1 THR    16     117.687  74.622  99.740  1.00  0.00    -0.393 OA
ATOM    180  HG1 THR    16     118.479  75.164  99.771  1.00  0.00     0.210 HD
ATOM    181  C   THR    16     115.992  72.841  96.911  1.00  0.00     0.243 C 
ATOM    182  O   THR    16     116.689  71.909  96.504  1.00  0.00    -0.271 OA
ATOM    183  N   GLY    17     115.063  73.437  96.167  1.00  0.00    -0.350 N 
ATOM    184  HN  GLY    17     114.498  74.180  96.552  1.00  0.00     0.163 HD
ATOM    185  CA  GLY    17     114.834  73.040  94.787  1.00  0.00     0.225 C 
ATOM    186  C   GLY    17     114.338  71.627  94.520  1.00  0.00     0.236 C 
ATOM    187  O   GLY    17     114.291  71.199  93.363  1.00  0.00    -0.272 OA
ATOM    188  N   ASP    18     113.966  70.888  95.560  1.00  0.00    -0.346 N 
ATOM    189  HN  ASP    18     114.027  71.230  96.508  1.00  0.00     0.163 HD
ATOM    190  CA  ASP    18     113.477  69.533  95.345  1.00  0.00     0.186 C 
ATOM    191  CB  ASP    18     113.305  68.797  96.682  1.00  0.00     0.147 C 
ATOM    192  CG  ASP    18     113.295  67.279  96.523  1.00  0.00     0.175 C 
ATOM    193  OD1 ASP    18     113.054  66.565  97.522  1.00  0.00    -0.648 OA
ATOM    194  OD2 ASP    18     113.536  66.797  95.395  1.00  0.00    -0.648 OA
ATOM    195  C   ASP    18     112.139  69.608  94.608  1.00  0.00     0.241 C 
ATOM    196  O   ASP    18     111.366  70.555  94.779  1.00  0.00    -0.271 OA
ATOM    197  N   LYS    19     111.883  68.621  93.761  1.00  0.00    -0.346 N 
ATOM    198  HN  LYS    19     112.571  67.903  93.583  1.00  0.00     0.163 HD
ATOM    199  CA  LYS    19     110.635  68.577  93.020  1.00  0.00     0.176 C 
ATOM    200  CB  LYS    19     110.863  68.947  91.545  1.00  0.00     0.035 C 
ATOM    201  CG  LYS    19     112.169  68.442  90.934  1.00  0.00     0.004 C 
ATOM    202  CD  LYS    19     112.255  68.822  89.447  1.00  0.00     0.027 C 
ATOM    203  CE  LYS    19     113.627  68.525  88.835  1.00  0.00     0.229 C 
ATOM    204  NZ  LYS    19     114.689  69.435  89.359  1.00  0.00    -0.079 N 
ATOM    205  HZ1 LYS    19     114.762  69.325  90.360  1.00  0.00     0.274 HD
ATOM    206  HZ2 LYS    19     115.573  69.203  88.929  1.00  0.00     0.274 HD
ATOM    207  HZ3 LYS    19     114.450  70.392  89.142  1.00  0.00     0.274 HD
ATOM    208  C   LYS    19     110.019  67.188  93.150  1.00  0.00     0.243 C 
ATOM    209  O   LYS    19     110.732  66.177  93.204  1.00  0.00    -0.271 OA
ATOM    210  N   PRO    20     108.681  67.123  93.228  1.00  0.00    -0.337 N 
ATOM    211  CD  PRO    20     107.725  68.238  93.104  1.00  0.00     0.127 C 
ATOM    212  CG  PRO    20     106.422  67.600  93.552  1.00  0.00     0.022 C 
ATOM    213  CB  PRO    20     106.535  66.204  93.030  1.00  0.00     0.037 C 
ATOM    214  CA  PRO    20     107.980  65.842  93.362  1.00  0.00     0.179 C 
ATOM    215  C   PRO    20     108.535  64.781  92.428  1.00  0.00     0.241 C 
ATOM    216  O   PRO    20     108.714  65.037  91.246  1.00  0.00    -0.271 OA
ATOM    229  N   MET    22     109.043  60.121  91.946  1.00  0.00    -0.346 N 
ATOM    230  HN  MET    22     109.214  60.250  90.959  1.00  0.00     0.163 HD
ATOM    231  CA  MET    22     108.787  58.773  92.432  1.00  0.00     0.177 C 
ATOM    232  CB  MET    22     108.031  57.950  91.382  1.00  0.00     0.045 C 
ATOM    233  CG  MET    22     106.680  58.505  90.975  1.00  0.00     0.076 C 
ATOM    234  SD  MET    22     105.533  58.606  92.361  1.00  0.00    -0.173 SA
ATOM    235  CE  MET    22     105.240  56.859  92.656  1.00  0.00     0.089 C 
ATOM    236  C   MET    22     110.098  58.064  92.743  1.00  0.00     0.241 C 
ATOM    237  O   MET    22     111.054  58.140  91.966  1.00  0.00    -0.271 OA
ATOM    238  N   GLN    23     110.131  57.373  93.878  1.00  0.00    -0.346 N 
ATOM    239  HN  GLN    23     109.314  57.358  94.471  1.00  0.00     0.163 HD
ATOM    240  CA  GLN    23     111.302  56.610  94.295  1.00  0.00     0.177 C 
ATOM    241  CB  GLN    23     111.814  57.109  95.650  1.00  0.00     0.044 C 
ATOM    242  CG  GLN    23     113.204  56.611  96.019  1.00  0.00     0.105 C 
ATOM    243  CD  GLN    23     113.759  57.245  97.297  1.00  0.00     0.215 C 
ATOM    244  OE1 GLN    23     113.694  58.467  97.491  1.00  0.00    -0.274 OA
ATOM    245  NE2 GLN    23     114.325  56.415  98.163  1.00  0.00    -0.370 N 
ATOM    246 1HE2 GLN    23     114.710  56.774  99.025  1.00  0.00     0.159 HD
ATOM    247 2HE2 GLN    23     114.369  55.427  97.957  1.00  0.00     0.159 HD
ATOM    248  C   GLN    23     110.815  55.172  94.412  1.00  0.00     0.241 C 
ATOM    249  O   GLN    23     109.811  54.907  95.074  1.00  0.00    -0.271 OA
ATOM    250  N   THR    24     111.507  54.248  93.756  1.00  0.00    -0.344 N 
ATOM    251  HN  THR    24     112.313  54.500  93.201  1.00  0.00     0.163 HD
ATOM    252  CA  THR    24     111.105  52.848  93.795  1.00  0.00     0.205 C 
ATOM    253  CB  THR    24     111.384  52.152  92.447  1.00  0.00     0.146 C 
ATOM    254  CG2 THR    24     110.669  50.794  92.384  1.00  0.00     0.042 C 
ATOM    255  OG1 THR    24     110.912  52.983  91.376  1.00  0.00    -0.393 OA
ATOM    256  HG1 THR    24     111.085  52.551  90.536  1.00  0.00     0.210 HD
ATOM    257  C   THR    24     111.841  52.100  94.897  1.00  0.00     0.243 C 
ATOM    258  O   THR    24     113.008  52.378  95.167  1.00  0.00    -0.271 OA
ATOM    259  N   TYR    25     111.147  51.167  95.544  1.00  0.00    -0.346 N 
ATOM    260  HN  TYR    25     110.186  50.987  95.293  1.00  0.00     0.163 HD
ATOM    261  CA  TYR    25     111.743  50.373  96.611  1.00  0.00     0.180 C 
ATOM    262  CB  TYR    25     111.206  50.772  97.998  1.00  0.00     0.073 C 
ATOM    263  CG  TYR    25     111.404  52.226  98.402  1.00  0.00    -0.056 A 
ATOM    264  CD1 TYR    25     110.388  53.168  98.220  1.00  0.00     0.010 A 
ATOM    265  CE1 TYR    25     110.561  54.505  98.598  1.00  0.00     0.037 A 
ATOM    266  CZ  TYR    25     111.760  54.912  99.167  1.00  0.00     0.065 A 
ATOM    267  OH  TYR    25     111.920  56.232  99.563  1.00  0.00    -0.361 OA
ATOM    268  HH  TYR    25     111.143  56.770  99.397  1.00  0.00     0.217 HD
ATOM    269  CE2 TYR    25     112.784  53.988  99.354  1.00  0.00     0.037 A 
ATOM    270  CD2 TYR    25     112.600  52.655  98.972  1.00  0.00     0.010 A 
ATOM    271  C   TYR    25     111.385  48.919  96.350  1.00  0.00     0.241 C 
ATOM    272  O   TYR    25     110.355  48.624  95.741  1.00  0.00    -0.271 OA
ATOM    273  N   GLU    26     112.233  48.014  96.818  1.00  0.00    -0.346 N 
ATOM    274  HN  GLU    26     113.082  48.309  97.279  1.00  0.00     0.163 HD
ATOM    275  CA  GLU    26     112.007  46.587  96.636  1.00  0.00     0.177 C 
ATOM    276  CB  GLU    26     113.149  45.995  95.801  1.00  0.00     0.045 C 
ATOM    277  CG  GLU    26     114.515  46.649  96.052  1.00  0.00     0.116 C 
ATOM    278  CD  GLU    26     115.554  46.293  94.996  1.00  0.00     0.172 C 
ATOM    279  OE1 GLU    26     115.353  46.658  93.817  1.00  0.00    -0.648 OA
ATOM    280  OE2 GLU    26     116.570  45.646  95.344  1.00  0.00    -0.648 OA
ATOM    281  C   GLU    26     111.897  45.923  98.007  1.00  0.00     0.241 C 
ATOM    282  O   GLU    26     112.612  46.286  98.938  1.00  0.00    -0.271 OA
ATOM    283  N   ILE    27     111.001  44.952  98.130  1.00  0.00    -0.346 N 
ATOM    284  HN  ILE    27     110.450  44.662  97.335  1.00  0.00     0.163 HD
ATOM    285  CA  ILE    27     110.766  44.293  99.413  1.00  0.00     0.180 C 
ATOM    286  CB  ILE    27     109.550  44.967 100.138  1.00  0.00     0.013 C 
ATOM    287  CG2 ILE    27     109.099  44.131 101.327  1.00  0.00     0.012 C 
ATOM    288  CG1 ILE    27     109.913  46.392 100.571  1.00  0.00     0.002 C 
ATOM    289  CD1 ILE    27     111.085  46.493 101.566  1.00  0.00     0.005 C 
ATOM    290  C   ILE    27     110.477  42.800  99.279  1.00  0.00     0.241 C 
ATOM    291  O   ILE    27     109.684  42.390  98.432  1.00  0.00    -0.271 OA
ATOM    292  N   ASP    28     111.110  41.986 100.119  1.00  0.00    -0.345 N 
ATOM    293  HN  ASP    28     111.763  42.346 100.800  1.00  0.00     0.163 HD
ATOM    294  CA  ASP    28     110.875  40.545 100.080  1.00  0.00     0.186 C 
ATOM    295  CB  ASP    28     112.019  39.784 100.765  1.00  0.00     0.147 C 
ATOM    296  CG  ASP    28     111.878  38.263 100.644  1.00  0.00     0.175 C 
ATOM    297  OD1 ASP    28     110.739  37.742 100.719  1.00  0.00    -0.648 OA
ATOM    298  OD2 ASP    28     112.915  37.582 100.484  1.00  0.00    -0.648 OA
ATOM    299  C   ASP    28     109.564  40.258 100.813  1.00  0.00     0.241 C 
ATOM    300  O   ASP    28     109.528  40.213 102.040  1.00  0.00    -0.271 OA
ATOM    301  N   LEU    29     108.489  40.060 100.060  1.00  0.00    -0.346 N 
ATOM    302  HN  LEU    29     108.550  40.173  99.058  1.00  0.00     0.163 HD
ATOM    303  CA  LEU    29     107.182  39.793 100.647  1.00  0.00     0.177 C 
ATOM    304  CB  LEU    29     106.150  39.547  99.545  1.00  0.00     0.038 C 
ATOM    305  CG  LEU    29     105.862  40.751  98.643  1.00  0.00    -0.020 C 
ATOM    306  CD1 LEU    29     104.632  40.454  97.787  1.00  0.00     0.009 C 
ATOM    307  CD2 LEU    29     105.641  41.999  99.497  1.00  0.00     0.009 C 
ATOM    308  C   LEU    29     107.109  38.660 101.669  1.00  0.00     0.241 C 
ATOM    309  O   LEU    29     106.093  38.510 102.349  1.00  0.00    -0.271 OA
ATOM    310  N   ASN    30     108.158  37.852 101.781  1.00  0.00    -0.346 N 
ATOM    311  HN  ASN    30     108.966  37.940 101.182  1.00  0.00     0.163 HD
ATOM    312  CA  ASN    30     108.143  36.775 102.772  1.00  0.00     0.185 C 
ATOM    313  CB  ASN    30     109.077  35.634 102.365  1.00  0.00     0.137 C 
ATOM    314  CG  ASN    30     108.337  34.483 101.711  1.00  0.00     0.217 C 
ATOM    315  OD1 ASN    30     107.976  34.548 100.537  1.00  0.00    -0.274 OA
ATOM    316  ND2 ASN    30     108.091  33.425 102.479  1.00  0.00    -0.370 N 
ATOM    317 1HD2 ASN    30     107.601  32.629 102.095  1.00  0.00     0.159 HD
ATOM    318 2HD2 ASN    30     108.395  33.422 103.442  1.00  0.00     0.159 HD
ATOM    319  C   ASN    30     108.564  37.345 104.126  1.00  0.00     0.241 C 
ATOM    320  O   ASN    30     108.241  36.791 105.180  1.00  0.00    -0.271 OA
ATOM    321  N   ASN    31     109.284  38.463 104.071  1.00  0.00    -0.346 N 
ATOM    322  HN  ASN    31     109.501  38.850 103.163  1.00  0.00     0.163 HD
ATOM    323  CA  ASN    31     109.762  39.172 105.260  1.00  0.00     0.185 C 
ATOM    324  CB  ASN    31     111.215  39.637 105.057  1.00  0.00     0.137 C 
ATOM    325  CG  ASN    31     112.191  38.486 104.880  1.00  0.00     0.217 C 
ATOM    326  OD1 ASN    31     113.327  38.688 104.454  1.00  0.00    -0.274 OA
ATOM    327  ND2 ASN    31     111.756  37.276 105.216  1.00  0.00    -0.370 N 
ATOM    328 1HD2 ASN    31     112.367  36.478 105.117  1.00  0.00     0.159 HD
ATOM    329 2HD2 ASN    31     110.817  37.160 105.569  1.00  0.00     0.159 HD
ATOM    330  C   ASN    31     108.877  40.406 105.486  1.00  0.00     0.241 C 
ATOM    331  O   ASN    31     109.371  41.475 105.855  1.00  0.00    -0.271 OA
ATOM    332  N   CYS    32     107.571  40.263 105.269  1.00  0.00    -0.345 N 
ATOM    333  HN  CYS    32     107.179  39.384 104.962  1.00  0.00     0.163 HD
ATOM    334  CA  CYS    32     106.676  41.400 105.425  1.00  0.00     0.185 C 
ATOM    335  CB  CYS    32     106.674  42.222 104.136  1.00  0.00     0.105 C 
ATOM    336  SG  CYS    32     105.933  43.855 104.290  1.00  0.00    -0.180 SA
ATOM    337  HG  CYS    32     106.123  44.214 103.024  1.00  0.00     0.101 HD
ATOM    338  C   CYS    32     105.248  41.017 105.769  1.00  0.00     0.241 C 
ATOM    339  O   CYS    32     104.872  39.842 105.736  1.00  0.00    -0.271 OA
ATOM    340  N   GLY    33     104.460  42.031 106.110  1.00  0.00    -0.351 N 
ATOM    341  HN  GLY    33     104.821  42.973 106.151  1.00  0.00     0.163 HD
ATOM    342  CA  GLY    33     103.066  41.816 106.432  1.00  0.00     0.225 C 
ATOM    343  C   GLY    33     102.275  41.747 105.141  1.00  0.00     0.238 C 
ATOM    344  O   GLY    33     102.836  41.951 104.060  1.00  0.00    -0.272 OA
ATOM    345  N   PRO    34     100.967  41.468 105.217  1.00  0.00    -0.337 N 
ATOM    346  CD  PRO    34     100.281  41.002 106.429  1.00  0.00     0.127 C 
ATOM    347  CG  PRO    34      99.442  39.897 105.877  1.00  0.00     0.022 C 
ATOM    348  CB  PRO    34      98.920  40.497 104.572  1.00  0.00     0.037 C 
ATOM    349  CA  PRO    34     100.076  41.363 104.057  1.00  0.00     0.179 C 
ATOM    350  C   PRO    34      99.556  42.694 103.545  1.00  0.00     0.241 C 
ATOM    351  O   PRO    34      99.151  42.795 102.390  1.00  0.00    -0.271 OA
ATOM    352  N   MET    35      99.564  43.711 104.401  1.00  0.00    -0.346 N 
ATOM    353  HN  MET    35      99.877  43.576 105.352  1.00  0.00     0.163 HD
ATOM    354  CA  MET    35      99.039  45.017 104.025  1.00  0.00     0.177 C 
ATOM    355  CB  MET    35      98.253  45.614 105.194  1.00  0.00     0.045 C 
ATOM    356  CG  MET    35      97.202  44.673 105.763  1.00  0.00     0.076 C 
ATOM    357  SD  MET    35      96.042  44.045 104.519  1.00  0.00    -0.173 SA
ATOM    358  CE  MET    35      94.846  43.170 105.591  1.00  0.00     0.089 C 
ATOM    359  C   MET    35     100.073  46.013 103.548  1.00  0.00     0.241 C 
ATOM    360  O   MET    35     101.252  45.921 103.875  1.00  0.00    -0.271 OA
ATOM    361  N   VAL    36      99.606  46.978 102.769  1.00  0.00    -0.346 N 
ATOM    362  HN  VAL    36      98.629  46.981 102.514  1.00  0.00     0.163 HD
ATOM    363  CA  VAL    36     100.465  48.009 102.214  1.00  0.00     0.180 C 
ATOM    364  CB  VAL    36      99.641  48.978 101.322  1.00  0.00     0.009 C 
ATOM    365  CG1 VAL    36     100.512  50.144 100.842  1.00  0.00     0.012 C 
ATOM    366  CG2 VAL    36      99.075  48.223 100.125  1.00  0.00     0.012 C 
ATOM    367  C   VAL    36     101.202  48.810 103.284  1.00  0.00     0.241 C 
ATOM    368  O   VAL    36     102.339  49.242 103.069  1.00  0.00    -0.271 OA
ATOM    369  N   LEU    37     100.559  49.014 104.432  1.00  0.00    -0.346 N 
ATOM    370  HN  LEU    37      99.623  48.659 104.569  1.00  0.00     0.163 HD
ATOM    371  CA  LEU    37     101.174  49.785 105.513  1.00  0.00     0.177 C 
ATOM    372  CB  LEU    37     100.169  50.062 106.636  1.00  0.00     0.038 C 
ATOM    373  CG  LEU    37     100.770  50.802 107.842  1.00  0.00    -0.020 C 
ATOM    374  CD1 LEU    37     101.237  52.208 107.437  1.00  0.00     0.009 C 
ATOM    375  CD2 LEU    37      99.733  50.867 108.948  1.00  0.00     0.009 C 
ATOM    376  C   LEU    37     102.367  49.043 106.080  1.00  0.00     0.241 C 
ATOM    377  O   LEU    37     103.346  49.658 106.509  1.00  0.00    -0.271 OA
ATOM    378  N   ASP    38     102.273  47.717 106.085  1.00  0.00    -0.346 N 
ATOM    379  HN  ASP    38     101.435  47.270 105.742  1.00  0.00     0.163 HD
ATOM    380  CA  ASP    38     103.350  46.875 106.577  1.00  0.00     0.186 C 
ATOM    381  CB  ASP    38     102.963  45.410 106.459  1.00  0.00     0.147 C 
ATOM    382  CG  ASP    38     101.773  45.055 107.318  1.00  0.00     0.175 C 
ATOM    383  OD1 ASP    38     101.789  45.413 108.516  1.00  0.00    -0.648 OA
ATOM    384  OD2 ASP    38     100.830  44.410 106.803  1.00  0.00    -0.648 OA
ATOM    385  C   ASP    38     104.596  47.144 105.746  1.00  0.00     0.241 C 
ATOM    386  O   ASP    38     105.706  47.268 106.279  1.00  0.00    -0.271 OA
ATOM    387  N   ALA    39     104.398  47.239 104.436  1.00  0.00    -0.346 N 
ATOM    388  HN  ALA    39     103.471  47.112 104.055  1.00  0.00     0.163 HD
ATOM    389  CA  ALA    39     105.490  47.505 103.517  1.00  0.00     0.172 C 
ATOM    390  CB  ALA    39     104.991  47.456 102.090  1.00  0.00     0.042 C 
ATOM    391  C   ALA    39     106.074  48.873 103.816  1.00  0.00     0.240 C 
ATOM    392  O   ALA    39     107.287  49.022 103.918  1.00  0.00    -0.271 OA
ATOM    393  N   LEU    40     105.206  49.870 103.962  1.00  0.00    -0.346 N 
ATOM    394  HN  LEU    40     104.217  49.699 103.850  1.00  0.00     0.163 HD
ATOM    395  CA  LEU    40     105.645  51.233 104.256  1.00  0.00     0.177 C 
ATOM    396  CB  LEU    40     104.431  52.142 104.415  1.00  0.00     0.038 C 
ATOM    397  CG  LEU    40     103.624  52.373 103.137  1.00  0.00    -0.020 C 
ATOM    398  CD1 LEU    40     102.196  52.740 103.485  1.00  0.00     0.009 C 
ATOM    399  CD2 LEU    40     104.292  53.469 102.312  1.00  0.00     0.009 C 
ATOM    400  C   LEU    40     106.480  51.251 105.527  1.00  0.00     0.241 C 
ATOM    401  O   LEU    40     107.571  51.812 105.549  1.00  0.00    -0.271 OA
ATOM    402  N   ILE    41     105.963  50.633 106.585  1.00  0.00    -0.346 N 
ATOM    403  HN  ILE    41     105.058  50.191 106.508  1.00  0.00     0.163 HD
ATOM    404  CA  ILE    41     106.665  50.557 107.864  1.00  0.00     0.180 C 
ATOM    405  CB  ILE    41     105.847  49.751 108.899  1.00  0.00     0.013 C 
ATOM    406  CG2 ILE    41     106.631  49.608 110.193  1.00  0.00     0.012 C 
ATOM    407  CG1 ILE    41     104.513  50.445 109.171  1.00  0.00     0.002 C 
ATOM    408  CD1 ILE    41     104.617  51.865 109.756  1.00  0.00     0.005 C 
ATOM    409  C   ILE    41     108.011  49.858 107.680  1.00  0.00     0.241 C 
ATOM    410  O   ILE    41     109.025  50.291 108.234  1.00  0.00    -0.271 OA
ATOM    411  N   LYS    42     108.008  48.775 106.902  1.00  0.00    -0.346 N 
ATOM    412  HN  LYS    42     107.140  48.425 106.522  1.00  0.00     0.163 HD
ATOM    413  CA  LYS    42     109.227  48.020 106.647  1.00  0.00     0.176 C 
ATOM    414  CB  LYS    42     108.921  46.748 105.857  1.00  0.00     0.035 C 
ATOM    415  CG  LYS    42     110.134  45.841 105.691  1.00  0.00     0.004 C 
ATOM    416  CD  LYS    42     109.785  44.533 104.979  1.00  0.00     0.027 C 
ATOM    417  CE  LYS    42     111.032  43.684 104.705  1.00  0.00     0.229 C 
ATOM    418  NZ  LYS    42     111.646  43.133 105.953  1.00  0.00    -0.079 N 
ATOM    419  HZ1 LYS    42     110.974  42.548 106.429  1.00  0.00     0.274 HD
ATOM    420  HZ2 LYS    42     112.460  42.585 105.714  1.00  0.00     0.274 HD
ATOM    421  HZ3 LYS    42     111.919  43.894 106.558  1.00  0.00     0.274 HD
ATOM    422  C   LYS    42     110.299  48.836 105.924  1.00  0.00     0.241 C 
ATOM    423  O   LYS    42     111.431  48.886 106.388  1.00  0.00    -0.271 OA
ATOM    424  N   ILE    43     109.967  49.479 104.805  1.00  0.00    -0.346 N 
ATOM    425  HN  ILE    43     109.049  49.405 104.390  1.00  0.00     0.163 HD
ATOM    426  CA  ILE    43     110.982  50.270 104.109  1.00  0.00     0.180 C 
ATOM    427  CB  ILE    43     110.528  50.737 102.674  1.00  0.00     0.013 C 
ATOM    428  CG2 ILE    43     109.425  49.853 102.145  1.00  0.00     0.012 C 
ATOM    429  CG1 ILE    43     110.090  52.197 102.687  1.00  0.00     0.002 C 
ATOM    430  CD1 ILE    43     111.198  53.163 102.321  1.00  0.00     0.005 C 
ATOM    431  C   ILE    43     111.339  51.477 104.972  1.00  0.00     0.241 C 
ATOM    432  O   ILE    43     112.385  52.100 104.791  1.00  0.00    -0.271 OA
ATOM    433  N   LYS    44     110.469  51.791 105.927  1.00  0.00    -0.346 N 
ATOM    434  HN  LYS    44     109.624  51.246 106.023  1.00  0.00     0.163 HD
ATOM    435  CA  LYS    44     110.684  52.907 106.850  1.00  0.00     0.176 C 
ATOM    436  CB  LYS    44     109.380  53.229 107.586  1.00  0.00     0.035 C 
ATOM    437  CG  LYS    44     109.512  54.184 108.767  1.00  0.00     0.004 C 
ATOM    438  CD  LYS    44     109.737  55.603 108.317  1.00  0.00     0.027 C 
ATOM    439  CE  LYS    44     109.715  56.562 109.501  1.00  0.00     0.229 C 
ATOM    440  NZ  LYS    44     109.912  57.976 109.055  1.00  0.00    -0.079 N 
ATOM    441  HZ1 LYS    44     109.172  58.234 108.417  1.00  0.00     0.274 HD
ATOM    442  HZ2 LYS    44     109.893  58.588 109.858  1.00  0.00     0.274 HD
ATOM    443  HZ3 LYS    44     110.804  58.061 108.588  1.00  0.00     0.274 HD
ATOM    444  C   LYS    44     111.762  52.544 107.871  1.00  0.00     0.241 C 
ATOM    445  O   LYS    44     112.657  53.338 108.153  1.00  0.00    -0.271 OA
ATOM    446  N   ASN    45     111.672  51.333 108.411  1.00  0.00    -0.346 N 
ATOM    447  HN  ASN    45     110.891  50.737 108.178  1.00  0.00     0.163 HD
ATOM    448  CA  ASN    45     112.624  50.866 109.411  1.00  0.00     0.185 C 
ATOM    449  CB  ASN    45     111.945  49.847 110.321  1.00  0.00     0.137 C 
ATOM    450  CG  ASN    45     110.759  50.428 111.038  1.00  0.00     0.217 C 
ATOM    451  OD1 ASN    45     110.713  51.627 111.305  1.00  0.00    -0.274 OA
ATOM    452  ND2 ASN    45     109.794  49.580 111.373  1.00  0.00    -0.370 N 
ATOM    453 1HD2 ASN    45     108.974  49.918 111.857  1.00  0.00     0.159 HD
ATOM    454 2HD2 ASN    45     109.884  48.601 111.143  1.00  0.00     0.159 HD
ATOM    455  C   ASN    45     113.944  50.268 108.906  1.00  0.00     0.241 C 
ATOM    456  O   ASN    45     114.929  50.224 109.645  1.00  0.00    -0.271 OA
ATOM    457  N   GLU    46     113.978  49.819 107.658  1.00  0.00    -0.346 N 
ATOM    458  HN  GLU    46     113.154  49.838 107.075  1.00  0.00     0.163 HD
ATOM    459  CA  GLU    46     115.189  49.203 107.125  1.00  0.00     0.177 C 
ATOM    460  CB  GLU    46     114.873  47.765 106.717  1.00  0.00     0.045 C 
ATOM    461  CG  GLU    46     114.374  46.922 107.873  1.00  0.00     0.116 C 
ATOM    462  CD  GLU    46     113.623  45.696 107.409  1.00  0.00     0.172 C 
ATOM    463  OE1 GLU    46     113.644  45.418 106.189  1.00  0.00    -0.648 OA
ATOM    464  OE2 GLU    46     113.014  45.011 108.262  1.00  0.00    -0.648 OA
ATOM    465  C   GLU    46     115.857  49.926 105.957  1.00  0.00     0.241 C 
ATOM    466  O   GLU    46     117.075  50.141 105.963  1.00  0.00    -0.271 OA
ATOM    467  N   ILE    47     115.066  50.306 104.962  1.00  0.00    -0.346 N 
ATOM    468  HN  ILE    47     114.076  50.109 104.989  1.00  0.00     0.163 HD
ATOM    469  CA  ILE    47     115.607  50.963 103.781  1.00  0.00     0.180 C 
ATOM    470  CB  ILE    47     114.722  50.700 102.566  1.00  0.00     0.013 C 
ATOM    471  CG2 ILE    47     115.491  51.023 101.295  1.00  0.00     0.012 C 
ATOM    472  CG1 ILE    47     114.221  49.251 102.593  1.00  0.00     0.002 C 
ATOM    473  CD1 ILE    47     115.311  48.174 102.744  1.00  0.00     0.005 C 
ATOM    474  C   ILE    47     115.801  52.467 103.882  1.00  0.00     0.241 C 
ATOM    475  O   ILE    47     116.894  52.975 103.645  1.00  0.00    -0.271 OA
ATOM    476  N   ASP    48     114.733  53.182 104.214  1.00  0.00    -0.345 N 
ATOM    477  HN  ASP    48     113.838  52.736 104.356  1.00  0.00     0.163 HD
ATOM    478  CA  ASP    48     114.804  54.633 104.313  1.00  0.00     0.186 C 
ATOM    479  CB  ASP    48     114.604  55.240 102.918  1.00  0.00     0.147 C 
ATOM    480  CG  ASP    48     114.403  56.752 102.946  1.00  0.00     0.175 C 
ATOM    481  OD1 ASP    48     114.759  57.409 103.950  1.00  0.00    -0.648 OA
ATOM    482  OD2 ASP    48     113.888  57.282 101.935  1.00  0.00    -0.648 OA
ATOM    483  C   ASP    48     113.785  55.190 105.302  1.00  0.00     0.241 C 
ATOM    484  O   ASP    48     112.585  55.270 105.008  1.00  0.00    -0.271 OA
ATOM    485  N   SER    49     114.275  55.573 106.476  1.00  0.00    -0.344 N 
ATOM    486  HN  SER    49     115.252  55.429 106.688  1.00  0.00     0.163 HD
ATOM    487  CA  SER    49     113.421  56.125 107.510  1.00  0.00     0.200 C 
ATOM    488  CB  SER    49     114.027  55.874 108.896  1.00  0.00     0.199 C 
ATOM    489  OG  SER    49     115.409  56.170 108.924  1.00  0.00    -0.398 OA
ATOM    490  HG  SER    49     115.754  56.002 109.804  1.00  0.00     0.209 HD
ATOM    491  C   SER    49     113.161  57.614 107.292  1.00  0.00     0.243 C 
ATOM    492  O   SER    49     112.764  58.328 108.216  1.00  0.00    -0.271 OA
ATOM    493  N   THR    50     113.396  58.086 106.071  1.00  0.00    -0.344 N 
ATOM    494  HN  THR    50     113.764  57.492 105.341  1.00  0.00     0.163 HD
ATOM    495  CA  THR    50     113.122  59.482 105.754  1.00  0.00     0.205 C 
ATOM    496  CB  THR    50     114.170  60.082 104.780  1.00  0.00     0.146 C 
ATOM    497  CG2 THR    50     113.683  61.401 104.189  1.00  0.00     0.042 C 
ATOM    498  OG1 THR    50     115.381  60.334 105.496  1.00  0.00    -0.393 OA
ATOM    499  HG1 THR    50     116.034  60.706 104.899  1.00  0.00     0.210 HD
ATOM    500  C   THR    50     111.749  59.495 105.106  1.00  0.00     0.243 C 
ATOM    501  O   THR    50     111.190  60.557 104.822  1.00  0.00    -0.271 OA
"""