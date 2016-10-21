import numpy

class Information():
    """A class for storing and accessing information about the elements of a
    pymolecule.Molecule object."""

    def __init__(self, parent_molecule_object):
        """Initializes the pymolecule.Information class.

            Args:
                parent_molecule_object -- The pymolecule.Molecule object
                    associated with this class.

            """

        self.__parent_molecule = parent_molecule_object

        self.__constants = {}
        #Removed HG from this list to avoid capturing gamma hydrogens
        self.__constants['element_names_with_two_letters'] = ['BR', 'CL', 'BI',
                                                              'AS', 'AG', 'LI',
                                                              'MG', 'RH', 'ZN',
                                                              'MN']

        #SHORTEN LENGTH OF BOND_LENGTH_DICT
        self.__constants['bond_length_dict'] = {
            'C-C': 1.53, 'N-N': 1.425, 'O-O': 1.469, 'S-S': 2.048,
            'C-H': 1.059, 'H-C': 1.059, 'C-N': 1.469, 'N-C': 1.469,
            'C-O': 1.413, 'O-C': 1.413, 'C-S': 1.819, 'S-C': 1.819,
            'N-H': 1.009, 'H-N': 1.009, 'N-O': 1.463, 'O-N': 1.463,
            'O-S': 1.577, 'S-O': 1.577, 'O-H': 0.967, 'H-O': 0.967,
            'S-H': 1.35, 'H-S': 1.35, 'S-N': 1.633, 'N-S': 1.633, 'C-F': 1.399,
            'F-C': 1.399, 'C-CL': 1.790, 'CL-C': 1.790, 'C-BR': 1.910,
            'BR-C': 1.910, 'C-I':2.162, 'I-C':2.162, 'S-BR': 2.321,
            'BR-S': 2.321, 'S-CL': 2.283, 'CL-S': 2.283, 'S-F': 1.640,
            'F-S': 1.640, 'S-I': 2.687, 'I-S': 2.687, 'P-BR': 2.366,
            'BR-P': 2.366, 'P-CL': 2.008, 'CL-P': 2.008, 'P-F': 1.495,
            'F-P': 1.495, 'P-I': 2.490, 'I-P': 2.490, 'P-C': 1.841,
            'C-P': 1.841, 'P-N': 1.730, 'N-P': 1.730, 'P-O': 1.662,
            'O-P': 1.662, 'P-S': 1.954, 'S-P': 1.954, 'N-BR': 1.843,
            'BR-N': 1.843, 'N-CL': 1.743, 'CL-N': 1.743, 'N-F': 1.406,
            'F-N': 1.406, 'N-I': 2.2, 'I-N': 2.2, 'SI-BR': 2.284,
            'BR-SI': 2.284, 'SI-CL': 2.072, 'CL-SI': 2.072, 'SI-F': 1.636,
            'F-SI': 1.636, 'SI-P': 2.264, 'P-SI': 2.264, 'SI-S': 2.145,
            'S-SI': 2.145, 'SI-SI': 2.359, 'SI-SI': 2.359, 'SI-C': 1.888,
            'C-SI': 1.888, 'SI-N': 1.743, 'N-SI': 1.743, 'SI-O': 1.631,
            'O-SI': 1.631, 'X-X': 1.53, 'X-C': 1.53, 'C-X': 1.53, 'X-H': 1.059,
            'H-X': 1.059, 'X-N': 1.469, 'N-X': 1.469, 'X-O': 1.413,
            'O-X': 1.413, 'X-S': 1.819, 'S-X': 1.819, 'X-F': 1.399,
            'F-X': 1.399, 'X-CL': 1.790, 'CL-X': 1.790, 'X-BR': 1.910,
            'BR-X': 1.910, 'X-I':2.162, 'I-X':2.162, 'SI-X': 1.888,
            'X-SI': 1.888, 'H-H':0.74}

        #INCLUDE ALL N AND C TERMINAL

        self.__constants['protein_residues'] = ["ALA", "ARG", "ASH", "ASN",
                                                "ASP", "CYM", "CYS", "CYX",
                                                "GLN", "GLH", "GLU", "GLY",
                                                "HID", "HIE", "HIP", "HIS",
                                                "ILE", "LEU", "LYN", "LYS",
                                                "MET", "PHE", "PRO", "SER",
                                                "THR", "TRP", "TYR", "VAL",
                                                "MSE", "TPO", "PTR", "SEP",
                                                'CALA', 'CARG', 'CASN', 'CASP',
                                                'CCYS', 'CCYX', 'CGLN', 'CGLU',
                                                'CGLY', 'CHID', 'CHIE', 'CHIP',
                                                'CHIS', 'CILE', 'CLEU', 'CLYS',
                                                'CMET', 'CPHE', 'CPRO', 'CSER',
                                                'CTHR', 'CTRP', 'CTYR', 'CVAL',
                                                'NALA', 'NARG', 'NASN', 'NASP',
                                                'NCYS', 'NCYX', 'NGLN', 'NGLU',
                                                'NGLY', 'NHID', 'NHIE', 'NHIP',
                                                'NHIS', 'NILE', 'NLEU', 'NLYS',
                                                'NMET', 'NPHE', 'NPRO', 'NSER',
                                                'NTHR', 'NTRP', 'NTYR', 'NVAL']

        self.__constants['dna_residues'] = ["A", "C", "G", "T", "DA", "DA3",
                                            "DA5", "DAN", "DC", "DC3", "DC4",
                                            "DC5", "DCN", "DG", "DG3", "DG5",
                                            "DGN", "DT", "DT3", "DT5", "DTN"]

        self.__constants['rna_residues'] = ["A", "C", "G", "U", "RA", "RA3",
                                            "RA5", "RAN", "RC", "RC3", "RC4",
                                            "RC5", "RCN", "RG", "RG3", "RG5",
                                            "RGN", "RU", "RU3", "RU5", "RUN"]

        self.__constants['mass_dict'] = {
            'H': 1.00794, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,
            'S': 32.065, 'P': 30.973762, 'NA': 22.9897, 'MG': 24.3050,
            'F': 18.9984032, 'CL': 35.453, 'K': 39.0983, 'CA': 40.078,
            'I': 126.90447, 'LI': 6.941, 'BE': 9.0122, 'B': 10.811,
            'AL': 26.9815, 'MN': 54.938, 'FE': 55.845, 'CO': 58.9332,
            'CU': 63.9332, 'ZN': 65.38, 'AS': 74.9216, 'BR': 79.904,
            'MO': 95.94, 'RH': 102.9055, 'AG': 107.8682, 'AU': 196.9655,
            'HG': 200.59, 'PB': 207.2, 'BI': 208.98040}

        self.__constants['vdw_dict'] = {
            'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.8,
            'S': 1.8, 'B': 2.0, 'LI': 1.82, 'NA': 2.27, 'MG': 1.73, 'AL': 2.00,
            'CL': 1.75, 'CA': 2.00, 'MN': 2.00, 'FE': 2.00, 'CO': 2.00,
            'CU': 1.40, 'ZN': 1.39, 'AS': 1.85, 'BR': 1.85, 'MO': 2.00,
            'RH': 2.00, 'AG': 1.72, 'AU': 1.66, 'PB': 2.02, 'BI': 2.00,
            'K': 2.75, 'I': 1.98
        }

        self.__constants['i8_fields'] = ['serial', 'resseq']
        self.__constants['f8_fields'] = ['x', 'y', 'z', 'occupancy',
                                        'tempfactor']
        self.__constants['max_number_of_bonds_permitted'] = {
            "C": 4, "N": 4, "O": 2, "H": 1, "F": 1, "Cl": 1, "BR": 1, "CL": 1,
            "I": 1, "P": 5, "S": 6
        }

        self.__filename = ""
        self.__remarks = []
        self.__atom_information = None
        self.__coordinates = None
        self.__coordinates_undo_point = None
        self.__bonds = None
        self.__hierarchy = {}
        self.__max_ring_size = 50

    #### Aliases ####
    # Gets
    def get_filename(self):
        return self.__filename

    def get_remarks(self):
        return self.__remarks

    def get_atom_information(self):
        return self.__atom_information

    def get_coordinates(self):
        return self.__coordinates

    def get_coordinates_undo_point(self):
        return self.__coordinates_undo_point

    def get_bonds(self):
        return self.__bonds

    def get_hierarchy(self):
        return self.__hierarchy

    def get_constants(self):
        return self.__constants

    # Sets
    def set_filename(self, filename):
        self.__filename = filename

    def set_remarks(self, remarks):
        self.__remarks = remarks

    def set_atom_information(self, atom_information):
        self.__atom_information = atom_information

    def set_coordinates(self, coordinates):
        self.__coordinates = coordinates

    def set_coordinates_undo_point(self, coordinates_undo_point):
        self.__coordinates_undo_point = coordinates_undo_point

    def set_bonds(self, bonds):
        self.__bonds = bonds

    def set_hierarchy(self, hierarchy):
        self.__hierarchy = hierarchy

    def belongs_to_protein(self, atom_index):
        """Checks if the atom is part of a protein. Taken primarily from Amber
        residue names.

            Args:
                atom_index -- An int, the index of the atom to consider.

            Returns:
                A boolean. True if part of protein, False if not.

            """

        # this function is retained for legacy reasons. past versions of
        # pymolecule had this functionality.

        if (self.__atom_information['resname_stripped'][atom_index]
            in self.__constants['protein_residues']):
            return True
        return False

    def belongs_to_dna(self, atom_index):
        """Checks if the atom is part of DNA.

            Args:
                atom_index -- An int, the index of the atom to consider.

            Returns:
                A boolean. True if part of dna, False if not.

            """

        # this function is retained for legacy reasons. past versions of
        # pymolecule had this functionality.

        if (self.__atom_information['resname_stripped'][atom_index]
            in self.__constants['dna_residues']):
            return True

        return False

    def belongs_to_rna(self, atom_index):
        """Checks if the atom is part of RNA.

            Args:
                atom_index -- An int, the index of the atom to consider.

            Returns:
                A boolean. True if part of rna, False if not.

            """

        # this function is retained for legacy reasons. past versions of
        # pymolecule had this functionality.

        if (self.__atom_information['resname_stripped'][atom_index]
            in self.__constants['rna_residues']):
            return True

        return False

    def assign_masses(self):
        """Assigns masses to the atoms of the pymolecule.Molecule object."""

        if not "mass" in self.__atom_information.dtype.names:
            # only assign if not been assigned previously
            masses = numpy.empty((
                len(self.__atom_information['element_stripped'])
            ))

            for i in range(len(self.__atom_information['element_stripped'])):
                element = self.__atom_information['element_stripped'][i]
                mass = self.__constants['mass_dict'][element]
                masses[i] = mass

            self.__atom_information = append_fields(self.__atom_information,
                                                    'mass', data = masses)

    def assign_elements_from_atom_names(self, selection = None):
        """Determines the elements of all atoms from the atom names. Note that
        this will overwrite any existing element assignments, including those
        explicitly specified in loaded files. Note that this doesn't populate
        elements_stripped.

            Args:
                selection -- An optional numpy.array containing the indices of
                    the atoms to consider when calculating the center of mass.
                    If ommitted, all atoms of the pymolecule.Molecule object
                    will be considered.

        """

        if selection is None:
            selection = self.__parent_molecule.select_all()

        if len(selection) == 0:
            return

        # get the atom names
        fix_element_names = numpy.core.defchararray.upper(
            self.__atom_information['name'][selection]
        )

        fix_element_names = numpy.core.defchararray.strip(fix_element_names)

        # first remove any numbers at the begining of these names
        fix_element_names = numpy.core.defchararray.lstrip(fix_element_names,
                                                           '0123456789')

        # remove any thing, letters or numbers, that follows a number,
        # including the number itself. so C2L becomes C, not CL.
        for num in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
            # I wish there was a more numpified way of doing this. :(
            tmp = numpy.core.defchararray.split(fix_element_names, num)
            fix_element_names = numpy.empty(len(fix_element_names),
                                            dtype = "S5")
            for i, item in enumerate(tmp):
                fix_element_names[i] = tmp[i][0]

        # take just first two letters of each item
        fix_element_names = numpy.array(fix_element_names, dtype = "|S2")

        # identify ones that are two-letter elements and one-letter elements
        cnsts = self.__constants
        one_tht_shf_b_2_lttrs = (
            fix_element_names == cnsts['element_names_with_two_letters'][0]
        )

        for other_two_letter in cnsts['element_names_with_two_letters'][1:]:
            one_tht_shf_b_2_lttrs = numpy.logical_or(
                one_tht_shf_b_2_lttrs,
                (fix_element_names == other_two_letter)
            )

        indices_of_two_letter_elements = numpy.nonzero(
            one_tht_shf_b_2_lttrs
        )[0]

        indices_of_one_letter_elements = numpy.nonzero(
            numpy.logical_not(one_tht_shf_b_2_lttrs)
        )[0]

        # get ones that are one-letter elements
        fix_element_names[indices_of_one_letter_elements] = (
            numpy.core.defchararray.rjust(numpy.array(
                fix_element_names[indices_of_one_letter_elements],
                dtype = "|S1"
            ), 2)
        )

        # they should be capitalized for consistency
        fix_element_names = numpy.core.defchararray.upper(fix_element_names)

        # now map missing element names back
        self.__atom_information['element'][selection] = fix_element_names

        # element_stripped also needs to be updated try:
        # self.__parent_molecule.information.get_atom_information()
        # ['element_stripped'][selection] =
        # numpy.core.defchararray.strip(fix_element_names) except: # so
        # element_stripped hasn't been defined yet
        #    self.__parent_molecule.information.get_atom_information() =
        #    append_fields(self.__parent_molecule.
        #    information.get_atom_information(), 'element_stripped',
        #    data = numpy.core.defchararray.strip(
        #    self.__parent_molecule.information.
        #    get_atom_information()['element']))

    def get_center_of_mass(self, selection = None):
        """Determines the center of mass.

            Args:
                selection -- An optional numpy.array containing the indices of
                    the atoms to consider when calculating the center of mass.
                    If ommitted, all atoms of the pymolecule.Molecule object
                    will be considered.

            Returns:
                A numpy.array containing to the x, y, and z coordinates of the
                    center of mass.

            """

        if selection is None: selection = self.__parent_molecule.select_all()

        # make sure the masses have been asigned
        self.assign_masses()

        # calculate the center of mass

        # multiply each coordinate by its mass
        center_of_mass = self.__coordinates[selection] * numpy.vstack((
            self.__atom_information['mass'][selection],
            self.__atom_information['mass'][selection],
            self.__atom_information['mass'][selection]
        )).T

        # now sum all that
        center_of_mass = numpy.sum(center_of_mass, 0)

        # now divide by the total mass
        center_of_mass = center_of_mass / self.get_total_mass(selection)

        return center_of_mass

    def get_geometric_center(self, selection = None):
        """Determines the geometric center.

            Args:
                selection -- An optional numpy.array containing the indices of
                    the atoms to consider when calculating the geometric
                    center. If ommitted, all atoms of the pymolecule.Molecule
                    object will be considered.

            Returns:
                A numpy.array containing to the x, y, and z coordinates of the
                    geometric center.

            """

        if selection is None:
            selection = self.__parent_molecule.select_all()

        return (numpy.sum(self.__coordinates[selection], 0) /
                self.get_total_number_of_atoms(selection))

    def get_total_mass(self, selection = None):
        """Calculates the total atomic mass.

            Args:
                selection -- An optional numpy.array containing the indices of
                    the atoms to consider when calculating the total mass. If
                    ommitted, all atoms of the pymolecule.Molecule object will
                    be considered.

            Returns:
                A double, the total mass.

        """

        if selection is None:
            selection = self.__parent_molecule.select_all()

        # assign masses if necessary
        self.assign_masses()

        # return total mass
        return numpy.sum(self.__atom_information['mass'][selection])

    def get_total_number_of_atoms(self, selection = None):
        """Counts the number of atoms.

            Args:
                selection -- An optional numpy.array containing the indices of
                    the atoms to count. If ommitted, all atoms of the
                    pymolecule.Molecule object will be considered.

            Returns:
                An int, the total number of atoms.

        """

        if selection is None:
            selection = self.__parent_molecule.select_all()

        if self.__coordinates is None:
            return 0
        else:
            return len(self.__coordinates[selection])

    def get_total_number_of_heavy_atoms(self):
        """Counts the number of heavy atoms (i.e., atoms that are not
        hydrogens).

            Args:
                selection -- An optional numpy.array containing the indices of
                    the atoms to count. If ommitted, all atoms of the
                    pymolecule.Molecule object will be considered.

            Returns:
                An int, the total number of heavy (non-hydrogen) atoms.

        """

        if self.__coordinates is None:
            return 0

        all_hydrogens = self.__parent_molecule.select_atoms({
            'element_stripped':'H'
        })

        return self.get_total_number_of_atoms() - len(all_hydrogens)

    def get_bounding_box(self, selection = None, padding = 0.0):
        """Calculates a box that bounds (encompasses) a set of atoms.

            Args:
                selection -- An optional numpy.array containing the indices of
                    the atoms to consider. If ommitted, all atoms of the
                    pymolecule.Molecule object will be considered.
                padding -- An optional float. The bounding box will extend this
                    many angstroms beyond the atoms being considered.

            Returns:
                A numpy array representing two 3D points, (min_x, min_y, min_z)
                    and (max_x, max_y, max_z), that bound the molecule.

            """

        if selection is None: selection = self.__parent_molecule.select_all()

        return numpy.vstack((numpy.min(self.__coordinates[selection], 0),
                             numpy.max(self.__coordinates[selection], 0)))

    def get_bounding_sphere(self, selection = None, padding = 0.0):
        """Calculates a sphere that bounds (encompasses) a set of atoms.

            Args:
                selection -- An optional numpy.array containing the indices of
                    the atoms to consider. If ommitted, all atoms of the
                    pymolecule.Molecule object will be considered.
                padding -- An optional float. The bounding sphere will extend
                    this many angstroms beyond the atoms being considered.

            Returns:
                A tuple containing two elements. The first is a numpy.array
                    representing a 3D point, the (x, y, z) center of the
                    sphere. The second is a float, the radius of the sphere.

            """

        if selection is None:
            selection = self.__parent_molecule.select_all()

        # get center
        center_of_selection = numpy.array([
            self.get_geometric_center(selection)
        ])

        # get distance to farthest point in selection
        return (center_of_selection[0],
                numpy.max(
                    cdist(center_of_selection,
                          self.__coordinates[selection])[0])
                )

    def define_molecule_chain_residue_spherical_boundaries(self):
        """Identifies spheres that bound (encompass) the entire molecule, the
        chains, and the residues. This information is stored in
        pymolecule.Molecule.information.hierarchy."""

        # first, check to see if it's already been defined
        if 'spheres' in self.__hierarchy.keys():
            return

        # set up the new structure
        hrchy = self.__hierarchy
        hrchy['spheres'] = {}
        hrchy['spheres']['molecule'] = {}
        hrchy['spheres']['chains'] = {}
        hrchy['spheres']['residues'] = {}

        # get all the chains and residues
        chains = self.__parent_molecule.selections_of_chains()
        residues = self.__parent_molecule.selections_of_residues()

        # do calcs for the whole molcules
        whole_mol_calc = self.get_bounding_sphere()
        hrchy['spheres']['molecule']['center'] = numpy.array(
            [whole_mol_calc[0]]
        )

        hrchy['spheres']['molecule']['radius'] = whole_mol_calc[1]

        # do calcs for the chains
        hrchy['spheres']['chains']['keys'] = numpy.array(
            # numpy string array e.g. ['a', 'b', 'c']
            hrchy['chains']['indices'].keys()
        )

        hrchy['spheres']['chains']['centers'] = numpy.empty(
            (len(hrchy['spheres']['chains']['keys']), 3)
        )

        hrchy['spheres']['chains']['radii'] = numpy.empty(
            len(hrchy['spheres']['chains']['keys'])
        )

        for index, chainid in enumerate(hrchy['spheres']['chains']['keys']):
            asphere = self.get_bounding_sphere(
                selection = hrchy['chains']['indices'][chainid]
            )

            hrchy['spheres']['chains']['centers'][index][0] = asphere[0][0]
            hrchy['spheres']['chains']['centers'][index][1] = asphere[0][1]
            hrchy['spheres']['chains']['centers'][index][2] = asphere[0][2]
            hrchy['spheres']['chains']['radii'][index] = asphere[1]

        # do calcs for the residues
        hrchy['spheres']['residues']['keys'] = numpy.array(
            hrchy['residues']['indices'].keys()
        )

        hrchy['spheres']['residues']['centers'] = numpy.empty(
            (len(hrchy['spheres']['residues']['keys']), 3)
        )

        hrchy['spheres']['residues']['radii'] = numpy.empty(
            len(hrchy['spheres']['residues']['keys'])
        )

        for index, resid in enumerate(hrchy['spheres']['residues']['keys']):
            asphere = self.get_bounding_sphere(
                selection = hrchy['residues']['indices'][resid]
            )
            hrchy['spheres']['residues']['centers'][index][0] = asphere[0][0]
            hrchy['spheres']['residues']['centers'][index][1] = asphere[0][1]
            hrchy['spheres']['residues']['centers'][index][2] = asphere[0][2]
            hrchy['spheres']['residues']['radii'][index] = asphere[1]

    def serial_reindex(self):
        """Reindexes the serial field of the atoms in the molecule, starting
        with 1."""

        for i in range(len(self.__atom_information['serial'])):
            self.__atom_information['serial'][i] = i + 1

    def resseq_reindex(self):
        """Reindexes the resseq field of the atoms in the molecule, starting
        with 1."""

        keys = numpy.core.defchararray.add(
            self.__atom_information['resname_stripped'], '-'
        )

        keys = numpy.core.defchararray.add(
            keys,
            numpy.array([str(t) for t in self.__atom_information['resseq']])
        )

        keys = numpy.core.defchararray.add(keys, '-')

        keys = numpy.core.defchararray.add(
            keys, self.__atom_information['chainid_stripped']
        )

        keys2 = numpy.insert(keys, 0, '')[:-1]
        index_of_change = numpy.nonzero(numpy.logical_not(keys == keys2))[0]
        index_of_change = numpy.append(index_of_change,
                                       len(self.__atom_information))

        count = 1
        for t in range(len(index_of_change[:-1])):
            start = index_of_change[t]
            end = index_of_change[t + 1]

            self.__atom_information['resseq'][numpy.arange(
                start, end, 1, dtype = 'int'
            )] = count

            count = count + 1
