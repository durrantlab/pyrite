import numpy
import sys
from scipy.spatial.distance import cdist

class Selections():
    """A class for selecting atoms"""

    ######## selections ########
    def __init__(self, parent_molecule_object):
        """Initializes the pymolecule.Selections class.

            Args:
                parent_molecule_object -- The pymolecule.Molecule object
                    associated with this class.

        """

        self.__parent_molecule = parent_molecule_object

    def select_atoms(self, selection_criteria):
        """Select a set of atoms based on user-specified criteria.

            Args:
                selection_criteria -- An dictionary, where the keys correspond
                    to keys in the
                    self.__parent_molecule.information.get_atom_information()
                    structured numpy array, and the values are lists of
                    acceptable matches. The selection is a logical "AND"
                    between dictionary entries, but "OR" within the value lists
                    themselves. For example: {'atom':['CA', 'O'], 'chain':'A',
                    'resname':'PRO'} would select all atoms with the names CA
                    or O that are located in the PRO residues of chain A.

            Returns:
                A numpy.array containing the indices of the atoms of the
                    selection.

            >>> _sample_selection().select_atoms(
            ...     {'element_stripped':'C'}
            ... )
            array([1, 2, 4, 5, 6])

        """

        try:
            # start assuming everything is selected
            selection = numpy.ones(
                len(self.__parent_molecule.get_atom_information()),
                dtype = bool
            )

            for key in selection_criteria.keys():

                vals = selection_criteria[key]

                # make sure the vals are in a list
                # if it's a single value, put it in a list
                if not type(vals) is list and not type(vals) is tuple:
                    vals = [vals]

                # make sure the vals are in the right format
                const = self.__parent_molecule.get_constants()
                if key in const['f8_fields']:
                    vals = [float(v) for v in vals]
                elif key in const['i8_fields']:
                    vals = [int(v) for v in vals]
                else:
                    vals = [v.strip() for v in vals]

                # "or" all the vals together
                # start assuming nothing is selected
                atm_inf = self.__parent_molecule.get_atom_information()
                subselection = numpy.zeros(
                    len(atm_inf), dtype = bool
                )

                for val in vals:
                    subselection = numpy.logical_or(
                        subselection,
                        (atm_inf[key] == val)
                    )

                # now "and" that with everything else
                selection = numpy.logical_and(selection, subselection)

            # now get the indices of the selection
            return numpy.nonzero(selection)[0]
        except:
            print "ERROR: Could not make the selection. Existing fields:"
            print "\t" + ", ".join(atm_inf.dtype.names)
            sys.exit(0)

    def select_atoms_in_bounding_box(self, bounding_box):
        """Selects all the atoms that are within a bounding box.

            Args:
                bounding_box -- A 2x3 numpy.array containing the minimum and
                    maximum points of the bounding box. Example:
                    numpy.array(
                        [[min_x, min_y, min_z], [max_x, max_y, max_z]]
                    ).

            Returns:
                A numpy.array containing the indices of the atoms that are
                    within the bounding box.

            >>> _sample_selection().select_atoms_in_bounding_box(
            ...     numpy.array([
            ...         [-2.1, 66.0, 75.6], 
            ...         [-1.5, 68.0, 79.0]
            ...     ])
            ... )
            array([3, 4])


        """

        min_pt = bounding_box[0]
        max_pt = bounding_box[1]
        coordinates = self.__parent_molecule.get_coordinates()
        sel1 = numpy.nonzero((coordinates[:, 0] > min_pt[0]))[0]
        sel2 = numpy.nonzero((coordinates[:, 0] < max_pt[0]))[0]
        sel3 = numpy.nonzero((coordinates[:, 1] > min_pt[1]))[0]
        sel4 = numpy.nonzero((coordinates[:, 1] < max_pt[1]))[0]
        sel5 = numpy.nonzero((coordinates[:, 2] > min_pt[2]))[0]
        sel6 = numpy.nonzero((coordinates[:, 2] < max_pt[2]))[0]
        sel = numpy.intersect1d(sel1, sel2)
        sel = numpy.intersect1d(sel, sel3)
        sel = numpy.intersect1d(sel, sel4)
        sel = numpy.intersect1d(sel, sel5)
        sel = numpy.intersect1d(sel, sel6)

        return sel

    def select_all_atoms_bound_to_selection(self, selection):
        """Selects all the atoms that are bound to a user-specified selection.

            Args:
                selection -- A numpy.array containing the indices of the
                    user-specified selection.

            Returns:
                A numpy.array containing the indices of the atoms that are
                bound to the user-specified selection. Note that this new
                selection does not necessarily include the indices of the
                original user-specified selection.

            >>> _sample_selection().select_all_atoms_bound_to_selection(
            ...     _sample_selection().select_atoms(
            ...         {'name_stripped':'CA'}
            ...     )    
            ... )
            array([0, 2, 4])

        """

        if self.__parent_molecule.information.get_bonds() is None:
            print ("You need to define the bonds to use" +
                   "select_all_atoms_bound_to_selection().")
            return

        bonds_to_consider = self.__parent_molecule.get_bonds()[selection]

        return numpy.unique(numpy.nonzero(bonds_to_consider)[1])

    def select_branch(self, root_atom_index, directionality_atom_index):
        """Identify an isolated "branch" of a molecular model. Assumes the
        atoms with indices root_atom_index and directionality_atom_index are
        bound to one another and that the branch starts at root_atom_index one
        and "points" in the direction of directionality_atom_index.

        Args:
            root_atom_index -- An int, the index of the first atom in the
                branch (the "root").
            directionality_atom_index -- An int, the index of the second atom
                in the branch, used to establish directionality

        Returns:
            A numpy array containing the indices of the atoms of the branch.

        """

        # note that this function is mostly retained for legacy reasons. the
        # old version of pymolecule had a branch-identification function.

        if self.__parent_molecule.get_bonds() is None:
            print ("To identify atoms in the same molecule as the atoms of " +
                   "a selection, you need to define the bonds.")
            return

        # Make sure atoms are neighboring
        if (not directionality_atom_index in
            self.select_all_atoms_bound_to_selection(
                numpy.array([root_atom_index]))
            ):

            print ("The root and directionality atoms, with indices " +
                   str(root_atom_index) + " and " +
                   str(directionality_atom_index) +
                   ", respectively, are not neighboring atoms.")

            return

        # first, set up the two indices need to manage the growing list of
        # connected atoms. current_index is the index in the list that you're
        # currently considering
        current_index = 1

        # create an "empty" array to store the indices of the connected atoms
        # can't know ahead of time what size, so let's use a python list #
        # -99999 *
        # numpy.ones(len(self.__parent_molecule.information.get_coordinates()),
        # dtype=int) # assume initially that all the atoms belong to this
        # molecule. this list will be shortened, possibly, later if that
        # assumption is incorrect.
        indices_of_this_branch = [root_atom_index, directionality_atom_index]

        while True:
            # get all the neighbors of the current atom
            try:
                current_atom_index = indices_of_this_branch[current_index]
            except:
                # this error because you've reached the end of the larger
                # molecule
                break

            neighbors_indices = (
                self.__parent_molecule.select_all_atoms_bound_to_selection(
                    numpy.array([current_atom_index])
                )
            )

            # get the ones in neighbors_indices that are not in
            # indices_of_this_molecule
            new_ones = numpy.setdiff1d(
                neighbors_indices, indices_of_this_branch
            )

            indices_of_this_branch.extend(new_ones)

            # prepare to look at the next atom in the list
            current_index = current_index + 1

        return numpy.array(indices_of_this_branch)

    def select_atoms_from_same_molecule(self, selection):
        """Selects all the atoms that belong to the same molecule as a
        user-defined selection, assuming that the pymolecule.Molecule object
        actually contains multiple physically distinct molecules that are not
        bound to each other via covalent bonds.

            Args:
                selection -- A numpy.array containing the indices of the
                    user-defined selection.

            Returns:
                A numpy.array containing the indices of the atoms belonging to
                    the same molecules as the atoms of the user-defined
                    selection.

            >>> _sample_selection().select_atoms_from_same_molecule(
            ...     _sample_selection().select_atoms(
            ...         {'name_stripped':'CA'}
            ...     )    
            ... )
            array([0, 1, 2, 3, 4, 5, 6])

        """

        # If your "Molecule" object actually contains several molecules, this
        # one selects all the atoms from any molecule containing any atom in
        # the selection note that bonds must be defined

        if self.__parent_molecule.get_bonds() is None:
            print ("To identify atoms in the same molecule as the atoms of " +
                   "a selection, you need to define the bonds.")
            return

        indices = []
        for index in selection:

            # first, set up the two indices need to manage the growing list of
            # connected atoms.
            # current_index is the index in the list that you're currently
            # considering
            current_index = 0

            # create an "empty" array to store the indices of the connected
            # atoms
            # can't know ahead of time what size, so let's use a python list #
            # -99999 *
            # numpy.ones(len(self.__parent_molecule.information.
            # get_coordinates()), dtype=int) # assume initially that all the
            # atoms belong to this molecule. this list will be shortened,
            # possibly, later if that assumption is incorrect.
            indices_of_this_molecule = [index]

            while True:
                # get all the neighbors of the current atom
                try:
                    current_atom_index = (
                        indices_of_this_molecule[current_index]
                    )
                except:
                    # this error because you've reached the end of the larger
                    # molecule
                    break

                prnt = self.__parent_molecule
                neighbors_indices = prnt.select_all_atoms_bound_to_selection(
                    numpy.array([current_atom_index])
                )

                # get the ones in neighbors_indices that are not in
                # indices_of_this_molecule
                new_ones = numpy.setdiff1d(neighbors_indices,
                                           indices_of_this_molecule)

                indices_of_this_molecule.extend(new_ones)

                # prepare to look at the next atom in the list
                current_index = current_index + 1

            # indices_of_this_molecule =
            # indices_of_this_molecule[:current_index-1] # so the list is
            # prunes down.
            indices.append(indices_of_this_molecule)

        # now merge and remove redundancies
        return numpy.unique(numpy.hstack(indices))

    def selections_of_constituent_molecules(self):
        """Identifies the indices of atoms belonging to separate molecules,
        assuming that the pymolecule.Molecule object actually contains multiple
        physically distinct molecules that are not bound to each other via
        covalent bonds.

            Returns:
                A python list of numpy.array objects containing the indices of
                    the atoms belonging to each molecule of the composite
                    pymolecule.Molecule object.

        """

        # If your pymolecule.Molecule object contains multiple molecules (e.g.,
        # several chains), this will return a list of selections corresponding
        # to the atoms of each molecule.

        atoms_not_yet_considered = self.select_all()
        selections = []

        while len(atoms_not_yet_considered) > 0:
            # add the atoms in the same molecule as the first atom in
            # atoms_not_yet_considered
            this_molecule_atoms = self.select_atoms_from_same_molecule(
                numpy.array([atoms_not_yet_considered[0]])
            )

            # now remove these from the atoms_not_yet_considered list
            atoms_not_yet_considered = numpy.setxor1d(
                this_molecule_atoms, atoms_not_yet_considered, True
            )

            # save the atoms of this molecule
            selections.append(this_molecule_atoms)

        return selections

    def select_atoms_near_other_selection(self, selection, cutoff):
        """Selects all atoms that are near the atoms of a user-defined
        selection.

            Args:
                selection -- A numpy.array containing the indices of the
                    user-defined selection.
                cutoff -- A float, the distance cutoff (in Angstroms).

            Returns:
                A numpy.array containing the indices of all atoms near the
                    user-defined selection, not including the atoms of the
                    user-defined selection themselves.

            >>> _sample_selection().select_atoms_near_other_selection(
            ...     _sample_selection().select_atoms(
            ...         {'name_stripped':'CA'}
            ...     ),
            ...     2.0
            ... )
            array([0, 2, 4])

        """

        # note that this does not return a selection that includes the input
        # selection. merge selections as required to get a selection that also
        # includes the input.

        invert_selection = self.invert_selection(selection)

        prnt = self.__parent_molecule

        selection_coors = prnt.get_coordinates()[selection]
        inversion_coors = prnt.get_coordinates()[invert_selection]

        indices_of_nearby = invert_selection[numpy.unique(
            numpy.nonzero(cdist(inversion_coors, selection_coors) < cutoff)[0]
        )]

        return indices_of_nearby

    def select_atoms_in_same_residue(self, selection):
        """Selects all atoms that are in the same residue as any of the atoms
        of a user-defined seleciton. Residues are considered unique if they
        have a unique combination of resname, resseq, and chainid fields.

            Args:
                selection -- A numpy.array containing the indices of the
                    user-defined selection.

            Returns:
                A numpy.array containing the indices of all atoms in the same
                    residue as any of the atoms of the user-defined selection.

            >>> _sample_selection().select_atoms_in_same_residue(
            ...     _sample_selection().select_atoms(
            ...         {'name_stripped':'CA'}
            ...     )
            ... )
            array([0, 1, 2, 3, 4, 5, 6])

        """

        # get string ids representing the residues of all atoms
        atm_inf = self.__parent_molecule.get_atom_information()

        keys = numpy.core.defchararray.add(atm_inf['resname_stripped'], '-')

        keys = numpy.core.defchararray.add(
            keys, numpy.array([str(t) for t in atm_inf['resseq']])
        )

        keys = numpy.core.defchararray.add(keys, '-')

        keys = numpy.core.defchararray.add(keys, atm_inf['chainid_stripped'])

        # get the unique keys of the selection
        unique_keys_of_selection = numpy.unique(keys[selection])

        # now get all the atoms of these selection keys

        # the below works, but is slow for large systems
        #residues = self.__parent_molecule.selections_of_residues()
        #new_selection = numpy.array([], dtype=int)
        #for key in unique_keys_of_selection:
        #    print key
        #    new_selection = numpy.append(new_selection, residues[key])

        # let's use this instead, faster for large systems.
        new_selection = numpy.array([], dtype = int)
        for key in unique_keys_of_selection:
            new_selection = numpy.append(
                new_selection, numpy.nonzero(keys == key)[0]
            )

        return new_selection

    def invert_selection(self, selection):
        """Inverts a user-defined selection (i.e., identifies all atoms that
        are not in the seleciton).

            Args:
                selection -- A numpy.array containing the indices of the
                    user-defined selection.

            Returns:
                A numpy.array containing the indices of all atoms that are not
                    in the user-defined seleciton.

            >>> _sample_selection().invert_selection(
            ...     _sample_selection().select_atoms(
            ...         {'name_stripped':'CA'}
            ...     )
            ... )
            array([0, 2, 3, 4, 5, 6])

        """

        # selection is a list of atom indices
        all_atoms = numpy.arange(
            0, len(self.__parent_molecule.get_atom_information()),
            1, dtype = int
        )

        remaining_indicies = numpy.delete(all_atoms, selection)
        return remaining_indicies

    def select_all(self):
        """Selects all the atoms in a pymolecule.Molecule object.

            Returns:
                A numpy.array containing the indices of all atoms in the
                    pymolecule.Molecule object.

            >>> _sample_selection().select_all()
            array([0, 1, 2, 3, 4, 5, 6])

        """

        return self.select_atoms({})

    def select_close_atoms_from_different_molecules(self, other_mol, cutoff,
                                                    pairwise_comparison = True,
                                                    terminate_early = False):
        """Effectively detects steric clashes between self and another
        pymolecule.Molecule.

            Args:
                other_mol -- A pymolecule.Molecule object of the other
                    molecule.
                cutoff -- A float, the user-defined distance cutoff in
                    Angstroms.
                pairwise_comparison -- An optional boolean, whether or not to
                    perform a simple pairwise distance comparison (if True) or
                    to use a more sophisitcated method (if False). True by
                    default.
                terminate_early = An optional boolean, whether or not to stop
                    looking for steric clashes once one is found. False by
                    default.

            Returns:
                A tuple containing two elements. The first is a numpy.array
                    containing the indices of all nearby atoms from this
                    pymolecule.Molecule object (self). The second is a
                    numpy.array containing the indices of all nearby atoms from
                    the other molecule.

        """

        if pairwise_comparison == True:

            dists = cdist(self.__parent_molecule.get_coordinates(),
                          other_mol.get_coordinates())
            close_ones = numpy.nonzero(dists < cutoff)
            close_ones_from_mol_parent_molecule = numpy.unique(close_ones[0])
            close_ones_from_mol_other_mol = numpy.unique(close_ones[1])

            return (close_ones_from_mol_parent_molecule,
                    close_ones_from_mol_other_mol)

        else:
            # so do the more complex hierarchical comparison
            # first, do some quick and easy checks
            self_coordinates = self.__parent_molecule.get_coordinates()
            other_coordinates = other_mol.get_coordinates()
            slf_hrachy = self.__parent_molecule.get_hierarchy()
            othr_hrachy = other_mol.get_hierarchy()

            margin = numpy.array([cutoff, cutoff, cutoff])
            self_min = numpy.min(self_coordinates, 0) - margin
            other_mol_max = numpy.max(other_coordinates, 0) + margin

            if self_min[0] > other_mol_max[0]:
                return (numpy.array([]), numpy.array([]))

            if self_min[1] > other_mol_max[1]:
                return (numpy.array([]), numpy.array([]))

            if self_min[2] > other_mol_max[2]:
                return (numpy.array([]), numpy.array([]))

            self_max = numpy.max(self_coordinates, 0) + margin
            other_mol_min = numpy.min(other_coordinates, 0) - margin

            if other_mol_min[0] > self_max[0]:
                return (numpy.array([]), numpy.array([]))

            if other_mol_min[1] > self_max[1]:
                return (numpy.array([]), numpy.array([]))

            if other_mol_min[2] > self_max[2]:
                return (numpy.array([]), numpy.array([]))

            # now assign spheres to the whole molecule, the chains, the
            # residues note that this won't recalculate the data if it's
            # already been calculated.
            prnt = self.__parent_molecule
            prnt.define_molecule_chain_residue_spherical_boundaries()
            other_mol.define_molecule_chain_residue_spherical_boundaries()

            # if the whole molecules are too far away, give up
            self_cent = slf_hrachy['spheres']['molecule']['center']
            self_rad = slf_hrachy['spheres']['molecule']['radius']
            other_cent = othr_hrachy['spheres']['molecule']['center']
            other_rad = othr_hrachy['spheres']['molecule']['radius']
            mol_dist = numpy.linalg.norm(self_cent - other_cent)

            if mol_dist > self_rad + other_rad + cutoff:
                # the molecules are too far away to clash
                return (numpy.array([]), numpy.array([]))

            # check the chains
            chain_distances = cdist(
                slf_hrachy['spheres']['chains']['centers'],
                other_mol.get_hierarchy()['spheres']['chains']['centers']
            )

            sum1_matrix = numpy.hstack(
                [
                    numpy.array(
                        [slf_hrachy['spheres']['chains']['radii']]
                    ).T for t in
                    range(len(othr_hrachy['spheres']['chains']['radii']))
                ]
            )

            sum2_matrix = numpy.vstack(
                [
                    numpy.array(
                        [othr_hrachy['spheres']['chains']['radii']]
                    ) for t in
                    range(len(slf_hrachy['spheres']['chains']['radii']))
                ]
            )

            sum_matrix = sum1_matrix + sum2_matrix + cutoff

            indicies_of_clashing_chains = numpy.nonzero(chain_distances <
                                                        sum_matrix)

            if len(indicies_of_clashing_chains[0]) == 0:
                # the chains don't clash, so no atoms can either
                return (numpy.array([]), numpy.array([]))

            # check the residues
            residue_distances = cdist(
                slf_hrachy['spheres']['residues']['centers'],
                other_mol.get_hierarchy()['spheres']['residues']['centers']
            )

            sum1_matrix = numpy.hstack(
                [
                    numpy.array(
                        [slf_hrachy['spheres']['residues']['radii']]
                    ).T for t in
                    range(len(othr_hrachy['spheres']['residues']['radii']))
                ]
            )

            sum2_matrix = numpy.vstack(
                [
                    numpy.array(
                        [othr_hrachy['spheres']['residues']['radii']]
                    ) for t in
                    range(len(slf_hrachy['spheres']['residues']['radii']))
                ]
            )

            sum_matrix = sum1_matrix + sum2_matrix + cutoff

            indicies_of_clashing_residues = numpy.nonzero(residue_distances <
                                                          sum_matrix)

            if len(indicies_of_clashing_residues[0]) == 0:
                # the residues don't clash, so no atoms can either
                return (numpy.array([]), numpy.array([]))

            # now time to check the atoms
            self_close_atom_indices = numpy.array([], dtype = int)
            other_close_atom_indices = numpy.array([], dtype = int)

            for i in range(len(indicies_of_clashing_residues[0])):
                self_res_index = indicies_of_clashing_residues[0][i]
                other_res_index = indicies_of_clashing_residues[1][i]

                self_res_name = (
                    slf_hrachy['spheres']['residues']['keys'][self_res_index]
                )

                other_res_name = (
                    othr_hrachy['spheres']['residues']['keys'][other_res_index]
                )

                self_res_indicies = (
                    slf_hrachy['residues']['indices'][self_res_name]
                )

                other_res_indicies = (
                    othr_hrachy['residues']['indices'][other_res_name]
                )

                self_coors = self_coordinates[self_res_indicies]
                other_coors = other_coordinates[other_res_indicies]

                some_self_indices, some_other_indices = numpy.nonzero(
                    cdist(self_coors, other_coors) < cutoff
                )

                if len(some_self_indices) != 0 or len(some_other_indices) != 0:
                    # so there are some
                    self_close_atom_indices = numpy.append(
                        self_close_atom_indices,
                        self_res_indicies[some_self_indices]
                    )

                    other_close_atom_indices = numpy.append(
                        other_close_atom_indices,
                        other_res_indicies[some_other_indices]
                    )

                    if terminate_early == True:
                        # so don't keep looking once you've found something
                        return (self_close_atom_indices,
                                other_close_atom_indices)

            # so nothing was found in the end
            return (numpy.unique(self_close_atom_indices),
                    numpy.unique(other_close_atom_indices))

    def get_molecule_from_selection(self, selection, serial_reindex = True,
                                    resseq_reindex = False):
        """Creates a pymolecule.Molecule from a user-defined atom selection.

            Args:
                selection -- A numpy.array containing the indices of the atoms
                    in the user-defined selection.
                serial_reindex -- An optional boolean, whether or not to
                    reindex the atom serial fields. Default is True.
                resseq_reindex -- An optional boolean, whether or not to
                    reindex the atom resseq fields. Default is False.

            Returns:
                A pymolecule.Molecule object containing the atoms of the
                    user-defined selection.

            >>> _sample_selection().get_molecule_from_selection(
            ...     _sample_selection().select_atoms(
            ...         {'name_stripped':'CA'}
            ...     ),
            ...     False, False
            ... ).get_atom_information()
            array([0, 2, 3, 4, 5, 6])

            """

        from pymolecule import Molecule
        new_mol = Molecule()
        new_mol.set_coordinates(
            self.__parent_molecule.get_coordinates()[selection]
        )

        # try to get the undo coordinates as well, though they may not have
        # been  set
        try:
            new_mol.information.set_coordinates_undo_point(
                self.__parent_molecule.get_coordinates_undo_point()[selection]
            )
        except:
            new_mol.information.set_coordinates_undo_point(None)

        new_mol.set_atom_information(
            self.__parent_molecule.get_atom_information()[selection]
        )

        if not self.__parent_molecule.get_bonds() is None:
            new_mol.set_bonds(self.__parent_molecule.get_bonds()[selection])
            new_mol.set_bonds(new_mol.get_bonds()[:, selection])
        else:
            new_mol.set_bonds(None)

        # note that hierarchy will have to be recalculated

        if serial_reindex == True: new_mol.information.serial_reindex()
        if resseq_reindex == True: new_mol.information.resseq_reindex()
        return new_mol

    def selections_of_chains(self):
        """Identifies the atom selections of each chain.

            Returns:
                A dictionary. The keys of the dictionary correspond to the
                    chainids, and the values are numpy.array objects containing
                    the indices of the associated chain atoms.

        """

        prnt = self.__parent_molecule

        if not 'chains' in prnt.get_hierarchy().keys():
            # so it hasn't already been calculated
            unique_chainids = numpy.unique(
                prnt.get_atom_information()['chainid_stripped']
            )

            prnt.get_hierarchy()['chains'] = {}
            prnt.get_hierarchy()['chains']['indices'] = {}
            for chainid in unique_chainids:
                prnt.get_hierarchy()['chains']['indices'][chainid] = (
                    prnt.select_atoms({'chainid_stripped': chainid})
                )

        return prnt.get_hierarchy()['chains']['indices']

    def selections_of_residues(self):
        """Identifies the atom selections of each residue.

            Returns:
                A dictionary. The keys of this dictionary correspond to the
                    unique resname-resseq-chainid residue identifiers, and the
                    values are numpy.array objects containing the indices of
                    the associated residue atoms.

            >>> _sample_selection().selections_of_residues()
            {'VAL-1-A': array([0, 1, 2, 3, 4, 5, 6])}

            """

        prnt = self.__parent_molecule
        atm_inf = prnt.get_atom_information()

        if not 'residues' in prnt.get_hierarchy().keys() :
            # so it hasn't already been calculated

            keys = numpy.core.defchararray.add(
                atm_inf['resname_stripped'], '-'
            )

            keys = numpy.core.defchararray.add(
                keys, numpy.array([str(t) for t in atm_inf['resseq']])
            )

            keys = numpy.core.defchararray.add(keys, '-')

            keys = numpy.core.defchararray.add(
                keys, atm_inf['chainid_stripped']
            )

            unique_resnames = numpy.unique(keys)

            prnt.get_hierarchy()['residues'] = {}

            prnt.get_hierarchy()['residues']['indices'] = {}

            for key in unique_resnames:
                resname, resseq, chainid = key.split('-')
                resseq = int(resseq)

                prnt.get_hierarchy()['residues']['indices'][key] = (
                    prnt.select_atoms({
                        'chainid_stripped': chainid,
                        'resname_stripped': resname,
                        'resseq': resseq})
                )

        return prnt.get_hierarchy()['residues']['indices']

    def in_same_ring(self, index1, index2):
        """Determines if two atoms in a Molecule are in the same ring.

        Args:
            index1 -- index of the first atom to see if it is in the ring.
            index2 -- index of the second atom to see if it is in the ring.

        Returns:
            A set of vertices of a path that has been traversed.

        """

        if index1 == index2: return True

        paths = []

        paths = self.__ring_recursive_walk(index1, index2, [], 0)

        if len(paths) == 0:
            print "No paths found between two indices"
            return False

        #Remove paths that do not start or end at the correct location
        for path in paths:
            if not (index2 in path and index1 in path): paths.remove(path)

        #Now need to find the intersection between each combination of paths
        for path1, path2 in itertools.combinations(paths, 2):
            intersection = set(path1).intersection(set(path2))

            #Only two elements are in index1 and index2
            if (len(intersection) == 2 and
                index1 in intersection and
                index2 in intersection):

                return True

        return False

    def __ring_recursive_walk(self, start, end, already_crossed, ringsize):
        """Recursive helper method to traverse a graph to search for a circular
        subgraph.

        Args:
            start -- graph vertex that is currently being traversed
            end -- vertex to be reached
            alreadyCrossed -- list of vertices already processed
            ringsize - counts number of vertices traversed

        Returns:
            A set of vertices of a path that has been traversed

        """

        paths = []

        ring_size += 1
        already_crossed.append(start)

        #Base case 1: Second point is reached
        if start == end:
            paths.append(already_crossed)
            return paths

        #Base case 2: Max ring size is reached
        if ring_size >= self.__max_ring_size():
            paths.append(already_crossed)
            return paths

        # Base case 3: No new neighbors
        # Get a list of all the atoms that atom:index is connected to that
        # haven't been previously evaluated
        sel = self.self.__parent_molecule.select_all_atoms_bound_to_selection
        neighbors = sel(numpy.array([start]))[:]

        for neighbor_index in neighbors:
            if neighbor_index in already_crossed:
                neighbors.remove(neighbor_index)

        if len(neighbors) == 0:
            paths.append(already_crossed)
            return paths

        for neighbor in neighbors:
            paths.extend(
                self.__ring_recursive_walk(neighbor, end,
                                           already_crossed[:], ring_size)
            )

        return paths


########## For unit testing ############

test_selection = None;

def _sample_selection():
    """A helper function that loads an example selection for testing."""

    global test_selection

    if test_selection is None:
        import Molecule as molecule
        import Selections as selections
        m = molecule.Molecule()
        m.load_pdb_into("./sample_structures/amino_acid_1.pdb", True, True, True)
        test_selection = selections.Selections(m)
    
    return test_selection

