from pymolecule import dumbpy as numpy
import os
import cPickle as pickle
import shutil

class FileIO():
    """A class for saving and loading molecular data into a pymolecule.Molecule
    object."""

    def __init__(self, parent_molecule_object):
        """Initializes the pymolecule.FileIO class.

            Args:
                parent_molecule_object -- The pymolecule.Molecule object
                    associated with this class.

        """

        self.__parent_molecule = parent_molecule_object

    def load_pym_into(self, filename):
        """Loads the molecular data contained in a pym file into the current
        pymolecule.Molecule object.

            Args:
                filename -- A string, the filename of the pym file.
        """

        if filename[-1:] != os.sep: filename = filename + os.sep

        # first, get the files that must exist
        self.__parent_molecule.set_atom_information(
            pickle.load(open(filename + 'atom_information', "rb"))
        )

        self.__parent_molecule.set_coordinates(
            numpy.load(filename + "coordinates.npz")['arr_0']
        )

        # now look for other possible files (optional output)
        prnt = self.__parent_molecule
        if os.path.exists(filename + 'remarks'):
            prnt.set_remarks(pickle.load(open(filename + 'remarks', "rb")))

        if os.path.exists(filename + 'hierarchy'):
            prnt.set_hierarchy(pickle.load(open(filename + 'hierarchy', "rb")))

        if os.path.exists(filename + 'filename'):
            prnt.set_filename(pickle.load(open(filename + 'filename', "rb")))
        
        if prnt.get_filename() == "":  # If still no filename, set it to the one used as a parameter.
            prnt.set_filename(filename)

        if os.path.exists(filename + "bonds.npz"):
            prnt.set_bonds(numpy.load(filename + "bonds.npz")['arr_0'])

        if os.path.exists(filename + "coordinates_undo_point.npz"):
            prnt.set_coordinates_undo_point(
                numpy.load(filename + "coordinates_undo_point.npz")['arr_0']
            )
        
    def load_pdbqt_into(self, filename, bonds_by_distance = False,
                      serial_reindex = True, resseq_reindex = False):
        """Loads the molecular data contained in a pdbqt file into the current
        pymolecule.Molecule object. Note that this implementation is
        incomplete. It doesn't save atomic charges, for example. The atom
        types are stored in the "element" and "element_stripped" columns.

            Args:
                filename -- A string, the filename of the pdbqt file.
                bonds_by_distance -- An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. False by
                    default, unlike for PDB.
                serial_reindex -- An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
                resseq_reindex -- An optional boolean, whether or not to
                    reindex the pdbqt resseq field. False by default.
        """

        self.__parent_molecule.set_filename(filename)

        # open/read the file
        # Treat PDBQT just like PDB, but merge last two columns.
        afile = open(filename, "r")
        self.load_pdbqt_into_using_file_object(afile, bonds_by_distance,
                                               serial_reindex, resseq_reindex)
        afile.close()

    def load_pdbqt_into_using_file_object(self, file_obj,
                                          bonds_by_distance = False,
                                          serial_reindex = True,
                                          resseq_reindex = False):
        """Loads molecular data from a python file object (pdbqt formatted)
        into the current pymolecule.Molecule object. Note that most users will
        want to use the load_pdb_into() function instead, which is identical
        except that it accepts a filename string instead of a python file
        object.

            Args:
                file_obj -- A python file object, containing pdb-formatted
                    data.
                bonds_by_distance -- An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. False by
                    default, unlike for PDB.
                serial_reindex -- An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
                resseq_reindex -- An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.

            """

        self.load_pdb_into_using_file_object(file_obj, bonds_by_distance,
                                             serial_reindex, resseq_reindex)

        # Now merge the last two columns.
        atom_inf = self.__parent_molecule.get_atom_information()
         
        atom_types = numpy.core.defchararray.add(
            atom_inf["element_stripped"], atom_inf["charge"]
        )

        atom_inf["element_stripped"] = numpy.defchararray_strip(atom_types)

        atom_inf["charge"] = "\n"

        atom_inf["element"] = numpy.core.defchararray.rjust(
            atom_inf["element_stripped"], 2
        )

        self.__parent_molecule.set_atom_information(atom_inf)

    def load_pdb_into(self, filename, bonds_by_distance = True,
                      serial_reindex = True, resseq_reindex = False):
        """Loads the molecular data contained in a pdb file into the current
        pymolecule.Molecule object.

            Args:
                filename -- A string, the filename of the pdb file.
                bonds_by_distance -- An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. True by
                    default.
                serial_reindex -- An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
                resseq_reindex -- An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
        """

        self.__parent_molecule.set_filename(filename)

        # open/read the file
        afile = open(filename, "r")
        self.load_pdb_into_using_file_object(afile, bonds_by_distance,
                                             serial_reindex, resseq_reindex)
        afile.close()

    def load_pdb_into_using_file_object(self, file_obj,
                                        bonds_by_distance = True,
                                        serial_reindex = True,
                                        resseq_reindex = False):
        """Loads molecular data from a python file object (pdb formatted) into
        the current pymolecule.Molecule object. Note that most users will want
        to use the load_pdb_into() function instead, which is identical except
        that it accepts a filename string instead of a python file object.

            Args:
                file_obj -- A python file object, containing pdb-formatted
                    data.
                bonds_by_distance -- An optional boolean, whether or not to
                    determine atomic bonds based on atom proximity. True by
                    default.
                serial_reindex -- An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
                resseq_reindex -- An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
            """

        # source_data = numpy.genfromtxt(file_obj,
        # dtype="S6,S5,S5,S4,S2,S4,S4,S8,S8,S8,S6,S6,S10,S2,S2",
        # names=['record_name', 'serial', 'name', 'resname', 'chainid',
        # 'resseq', 'empty', 'x', 'y', 'z', 'occupancy', 'tempfactor',
        # 'empty2', 'element', 'charge'], delimiter=[6, 5, 5, 4, 2, 4, 4, 8, 8,
        # 8, 6, 6, 10, 2, 2])

        source_data = numpy.genfromtxt(
            file_obj,
            dtype = "S6,S5,S5,S5,S1,S4,S4,S8,S8,S8,S6,S6,S10,S2,S3",
            names = ['record_name', 'serial', 'name', 'resname', 'chainid',
                   'resseq', 'empty', 'x', 'y', 'z', 'occupancy',
                   'tempfactor', 'empty2', 'element', 'charge'],
            delimiter = [6, 5, 5, 5, 1, 4, 4, 8, 8, 8, 6, 6, 10, 2, 3]
        )

        # get the remarks, if any. good to hold on to this because some of my
        # programs might retain info via remarks
        remark_indices = numpy.nonzero(
            source_data['record_name'] == "REMARK"
        )[0]

        remarks = []
        for index in remark_indices:
            astr = ""
            for name in source_data.dtype.names[1:]:
                astr = astr + source_data[name][index]
            remarks.append(astr.rstrip())
        self.__parent_molecule.set_remarks(remarks)

        # in case the pdb file has only one line
        if source_data.ndim == 0:
            source_data = source_data.reshape(1, -1)

        # get the ones that are ATOM or HETATOM in the record_name
        or_matrix = numpy.logical_or((source_data['record_name'] == "ATOM  "),
                                     (source_data['record_name'] == "HETATM"))
        indices_of_atom_or_hetatom = numpy.nonzero(or_matrix)[0]
        self.__parent_molecule.set_atom_information(
            source_data[indices_of_atom_or_hetatom]
        )

        # now, some of the data needs to change types
        # first, fields that should be numbers cannot be empty strings
        atom_inf = self.__parent_molecule.get_atom_information()
        for field in (self.__parent_molecule.get_constants()['i8_fields'] +
                      self.__parent_molecule.get_constants()['f8_fields']):
            check_fields = atom_inf[field]
            check_fields = numpy.defchararray_strip(check_fields)
            indices_of_empty = numpy.nonzero(check_fields == '')[0]
            atom_inf[field][indices_of_empty] = '0'

        # now actually change the type
        old_types = atom_inf.dtype
        descr = old_types.descr

        for field in self.__parent_molecule.get_constants()['i8_fields']:
            index = atom_inf.dtype.names.index(field)
            descr[index] = (descr[index][0], 'i8')
        for field in self.__parent_molecule.get_constants()['f8_fields']:
            index = atom_inf.dtype.names.index(field)
            descr[index] = (descr[index][0], 'f8')
        new_types = numpy.dtype(descr)
        self.__parent_molecule.set_atom_information(atom_inf.astype(new_types))

        # remove some of the fields that just contain empty data
        atom_inf = self.__parent_molecule.get_atom_information()
        self.__parent_molecule.set_atom_information(
            self.__parent_molecule.numpy_structured_array_remove_field(
                atom_inf, ['empty', 'empty2']
            )
        )

        # the coordinates need to be placed in their own special numpy array to
        # facilitate later manipulation
        atom_inf = self.__parent_molecule.get_atom_information()
        self.__parent_molecule.set_coordinates(
            numpy.vstack([atom_inf['x'], atom_inf['y'], atom_inf['z']]).T
        )

        # now remove the coordinates from the atom_information object to save
        # memory
        self.__parent_molecule.set_atom_information(
            self.__parent_molecule.numpy_structured_array_remove_field(
                atom_inf, ['x', 'y', 'z']
            )
        )

        # now determine element from atom name for those entries where it's not
        # given note that the
        # molecule.information.assign_elements_from_atom_names function can be
        # used to overwrite this and assign elements based on the atom name
        # only.
        indicies_where_element_is_not_defined = numpy.nonzero(
            numpy.defchararray_strip(atom_inf['element']) == ''
        )[0]

        self.__parent_molecule.assign_elements_from_atom_names(
            indicies_where_element_is_not_defined
        )

        # string values in
        # self.__parent_molecule.information.get_atom_information() should also
        # be provided in stripped format for easier comparison
        fields_to_strip = ['name', 'resname', 'chainid', 'element']
        for f in fields_to_strip:
            self.__parent_molecule.set_atom_information(
                numpy.append_fields(
                    self.__parent_molecule.get_atom_information().copy(),
                    f + '_stripped',
                    data = numpy.defchararray_strip(atom_inf[f])
                )
            )

        # now, if there's conect data, load it. this part of the code is not
        # that "numpyic"
        conect_indices = numpy.nonzero(
            source_data['record_name'] == "CONECT"
        )[0]

        if len(conect_indices) > 0:

            self.__parent_molecule.set_bonds(numpy.zeros((len(atom_inf),
                                                          len(atom_inf))))

            # build serial to index mapping
            serial_to_index = {}
            # is there a faster way?
            for index, inf in enumerate(atom_inf['serial']):
                serial_to_index[inf] = index

            # get the connect data
            for index in conect_indices:
                astr = ""
                for name in source_data.dtype.names[1:]:
                    astr = astr + source_data[name][index]
                astr = astr.rstrip()

                indices = []
                for i in xrange(0, len(astr), 5):
                    indices.append(serial_to_index[int(astr[i:i + 5])])

                for partner_index in indices[1:]:
                    self.__parent_molecule.add_bond(indices[0], partner_index)

        # else: # create empty bond array
        #    self.__parent_molecule.information.get_bonds() =
        #    numpy.zeros((len(self.__parent_molecule.information.
        #        get_atom_information()),
        #    len(self.__parent_molecule.information.get_atom_information())))

        if bonds_by_distance == True:
            self.__parent_molecule.create_bonds_by_distance(True)

        if serial_reindex == True:
            self.__parent_molecule.serial_reindex()

        if resseq_reindex == True:
            self.__parent_molecule.resseq_reindex()

    def save_pym(self, filename, save_bonds = False, save_filename = False,
                 save_remarks = False, save_hierarchy = False,
                 save_coordinates_undo_point = False):
        """Saves the molecular data contained in a pymolecule.Molecule object
        to a pym file.

            Args:
                filename -- An string, the filename to use for saving. (Note
                    that this is actually a directory, not a file.)
                save_bonds -- An optional boolean, whether or not to save
                    information about atomic bonds. False by default.
                save_filename -- An optional boolean, whether or not to save
                    the original (pdb) filename. False by default.
                save_remarks -- An optional boolean, whether or not to save
                    remarks associated with the molecule. False by default.
                save_hierarchy -- An optional boolean, whether or not to save
                    information about spheres the bound (encompass) the whole
                    molecule, the chains, and the residues. False by default.
                save_coordinates_undo_point -- An optional boolean, whether or
                    not to save the last coordinate undo point. False by
                    default.
        """

        # Why not just pickle self.parent.information? Because it's a huge
        # file, can't selectively not save bonds, for example, and numpy.save
        # is faster than cPickle protocol 2 on numpy arrays

        # if the directory already exists, first delete it
        if os.path.exists(filename):
            try:
                shutil.rmtree(filename)
            except:
                pass

            # it could be a file, not a directory
            try:
                os.remove(filename)
            except:
                pass

        # filename is actually a directory, so append separator if needed
        if filename[-1:] != os.sep:
            filename = filename + os.sep

        # make directory
        os.mkdir(filename)

        # save components

        # python objects must be pickled
        if save_hierarchy == True:
            # note this is a combo of python objects and numpy arrays, so must
            # be pickled.
            pickle.dump(self.__parent_molecule.get_hierarchy(),
                        open(filename + 'hierarchy', 'wb'), -1)

        if save_remarks == True:
            # using the latest protocol
            pickle.dump(self.__parent_molecule.get_remarks(),
                        open(filename + 'remarks', 'wb'), -1)

        if save_filename == True:
            pickle.dump(self.__parent_molecule.get_filename(),
                        open(filename + 'filename', 'wb'), -1)

        # unfortunately, the speedy numpy.save doesn't work on masked arrays
        # masked arrays have a dump method, but it just uses cPickle so we're
        # just going to cPickle masked arrays. Could be so much faster if numpy
        # were up to speed... :( not clear that numpy.ma.dump accepts protocol
        # parameter, so let's just use cPickle directly
        pickle.dump(self.__parent_molecule.get_atom_information(),
                    open(filename + 'atom_information', 'wb'), -1)

        # fortunately, coordinates and bonds are regular numpy arrays they can
        # be saved with numpy's speedy numpy.save function note that I'm
        # compressing them here. benchmarking suggests this takes longer to
        # save, but is much faster to load. so I'm prioritizing load times over
        # save times note also that numpy.savez can save multiple arrays to a
        # single file, probably speeding up load.

        numpy.savez(filename + "coordinates.npz",
                    self.__parent_molecule.get_coordinates())

        if save_bonds == True:
            numpy.savez(filename + "bonds.npz",
                        self.__parent_molecule.get_bonds())

        if save_coordinates_undo_point == True:
            numpy.savez(filename + "coordinates_undo_point.npz",
                        self.__parent_molecule.get_coordinates_undo_point())

    def save_pdb(self, filename = "", serial_reindex = True,
                 resseq_reindex = False, return_text = False):
        """Saves the molecular data contained in a pymolecule.Molecule object
        to a pdb file.

            Args:
                filename -- An string, the filename to use for saving.
                serial_reindex -- An optional boolean, whether or not to
                    reindex the pdb serial field. True by default.
                resseq_reindex -- An optional boolean, whether or not to
                    reindex the pdb resseq field. False by default.
                return_text -- An optional boolean, whether or not to return
                    text instead of writing to a file. If True, the filename
                    variable is ignored.

            Returns:
                If return_text is True, a PDB-formatted string. Otherwise,
                returns nothing.
        """

        if len(self.__parent_molecule.get_atom_information()) > 0:
            # so the pdb is not empty (if it is empty, don't save)

            if serial_reindex == True: self.__parent_molecule.serial_reindex()
            if resseq_reindex == True: self.__parent_molecule.resseq_reindex()

            if return_text == False:
                afile = open(filename, "w")
            else:
                return_string = ""

            # print out remarks
            for line in self.__parent_molecule.get_remarks():
                remark = "REMARK" + line + "\n"

                if return_text == False:
                    afile.write(remark)
                else:
                    return_string = return_string + remark

            # print out coordinates
            atom_information = self.__parent_molecule.get_atom_information()
            coordinates = self.__parent_molecule.get_coordinates()

            printout = numpy.core.defchararray.add(
                atom_information['record_name'],
                numpy.core.defchararray.rjust(
                    atom_information['serial'].astype('|S5'), 5
                )
            )

            printout = numpy.core.defchararray.add(printout,
                                                   atom_information['name'])

            printout = numpy.core.defchararray.add(printout,
                                                   atom_information['resname'])

            printout = numpy.core.defchararray.add(printout,
                                                   atom_information['chainid'])

            printout = numpy.core.defchararray.add(
                printout, numpy.core.defchararray.rjust(
                    atom_information['resseq'].astype('|S4'), 4
                )
            )

            printout = numpy.core.defchararray.add(printout, '    ')

            printout = numpy.core.defchararray.add(
                printout, numpy.core.defchararray.rjust(
                    numpy.array(["%.3f" % t for t in coordinates[:, 0]]), 8
                )
            )

            printout = numpy.core.defchararray.add(
                printout, numpy.core.defchararray.rjust(
                    numpy.array(["%.3f" % t for t in coordinates[:, 1]]), 8
                )
            )

            printout = numpy.core.defchararray.add(
                printout, numpy.core.defchararray.rjust(
                    numpy.array(["%.3f" % t for t in coordinates[:, 2]]), 8
                )
            )

            printout = numpy.core.defchararray.add(
                printout, numpy.core.defchararray.rjust(
                    numpy.array(["%.2f" % t
                                 for t in atom_information['occupancy']]),
                    6
                )
            )

            printout = numpy.core.defchararray.add(
                printout, numpy.core.defchararray.rjust(
                    numpy.array(["%.2f" % t
                                 for t in atom_information['tempfactor']]),
                    6
                )
            )

            printout = numpy.core.defchararray.add(printout, '          ')

            printout = numpy.core.defchararray.add(
                printout, atom_information['element']
            )

            printout = numpy.core.defchararray.add(
                printout, atom_information['charge']
            )

            if return_text == False:
                if printout[0][-1:] == "\n":
                    afile.write("".join(printout) + "\n")
                else:
                    afile.write("\n".join(printout) + "\n")
            else:
                if printout[0][-1:] == "\n":
                    return_string += "".join(printout) + "\n"
                else:
                    return_string += "\n".join(printout) + "\n"

            # print out connect
            prnt = self.__parent_molecule
            atm_inf = prnt.get_atom_information()
            sel_atms_bnd_to_sel = prnt.select_all_atoms_bound_to_selection
            if not prnt.get_bonds() is None:
                for indx in range(len(prnt.get_bonds())):
                    indices_of_bond_partners = sel_atms_bnd_to_sel(
                        numpy.array([indx])
                    )

                    if len(indices_of_bond_partners) > 0:

                        if return_text == False:
                            afile.write(
                                "CONECT" +
                                str(atm_inf["serial"][indx]).rjust(5) +
                                "".join([str(atm_inf["serial"][t]).rjust(5)
                                         for t in indices_of_bond_partners]) +
                                "\n"
                            )
                        else:
                            return_string += (
                                "CONECT" +
                                str(atm_inf["serial"][indx]).rjust(5) +
                                "".join([str(atm_inf["serial"][t]).rjust(5)
                                         for t in indices_of_bond_partners]) +
                                "\n"
                            )

            if return_text == False:
                afile.close()
            else:
                return return_string

        else:
            print ("ERROR: Cannot save a Molecule with no atoms " +
                   "(file name \"" + filename + "\")")
