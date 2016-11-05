from pymolecule import dumbpy as numpy

class Manipulation():
    """A class for translating and rotating the atomic coordinates of a
    pymolecule.Molecule object"""

    def __init__(self, parent_molecule_object):
        """Initializes the pymolecule.Manipulation class.

            Args:
                parent_molecule_object -- The pymolecule.Molecule object
                    associated with this class.

            """

        self.__parent_molecule = parent_molecule_object

    def set_coordinate_undo_point(self):
        """Sets ("saves") the undo point of the atom coordinates. Any
        subsequent manipulations of atomic coordinates can be "undone" by
        reseting to this configuration via the coordinate_undo function."""

        self.__parent_molecule.set_coordinates_undo_point(
            self.__parent_molecule.get_coordinates().copy()
        )

    def coordinate_undo(self):
        """Resets the coordinates of all atoms to those saved using the
        set_coordinate_undo_point function."""

        self.__parent_molecule.set_coordinates(
            self.__parent_molecule.get_coordinates_undo_point().copy()
        )

    def set_atom_location(self, atom_index, new_location):
        """Translates the entire molecular model (without rotating) so that the
        atom with the specified index is located at the specified coordinate.

            Args:
                atom_index -- An int, the index of the target atom.
                new_location -- A numpy.array specifying the new (x, y, z)
                    coordinate of the specified atom.

            Returns:
                A numpy.array specifying the (delta_x, delta_y, delta_z) vector
                by which the pmolecule.Molecule was translated.

            """

        if new_location.shape == (3,):
            new_location = numpy.array([new_location])

        currentloc = self.__parent_molecule.get_coordinates()[atom_index]
        delta = new_location - currentloc

        self.translate_molecule(delta)

        return delta

    def translate_molecule(self, delta):
        """Translate all the atoms of the molecular model by a specified
        vector.

        Args:
            delta -- A numpy.array (delta_x, delta_y, delta_z) specifying the
            amount to move each atom along the x, y, and z coordinates.

        """

        if delta.shape == (3,):
            delta = numpy.array([delta])

        self.__parent_molecule.set_coordinates(
            self.__parent_molecule.get_coordinates() + delta
        )

        if 'spheres' in self.__parent_molecule.get_hierarchy().keys():
            gt_hrchy_sph = self.__parent_molecule.get_hierarchy()['spheres']

            # so update location of hierarchical elements
            gt_hrchy_sph['molecule']['center'] = (
                gt_hrchy_sph['molecule']['center'] + delta
            )

            gt_hrchy_sph['chains']['centers'] = (
                gt_hrchy_sph['chains']['centers'] + delta
            )

            gt_hrchy_sph['residues']['centers'] = (
                gt_hrchy_sph['residues']['centers'] + delta
            )

    def rotate_molecule_around_a_line_between_points(self, line_point1,
                                                     line_point2, rotate):
        """Rotate the molecular model about a line segment. The end points of
        the line segment are explicitly specified coordinates.

            Args:
                line_point1 -- A numpy.array (x, y, z) corresponding to one end
                    of the line segment.
                line_point2 -- A numpy.array (x, y, z) corresponding to the
                    other end of the line segment.
                rotate -- A float, the angle of rotation, in radians.

            """

        if line_point1.shape == (1, 3):
            line_point1 = line_point1[0]

        if line_point2.shape == (1, 3):
            line_point2 = line_point2[0]

        a = line_point1[0]
        b = line_point1[1]
        c = line_point1[2]
        #d = line_point2[0]
        #e = line_point2[1]
        #f = line_point2[2]

        delta = line_point2 - line_point1

        u = delta[0] #d - a
        v = delta[1] #e - b
        w = delta[2] #f - c

        v_2_plus_w_2 = numpy.power(v, 2) + numpy.power(w, 2)
        u_2_plus_w_2 = numpy.power(u, 2) + numpy.power(w, 2)
        u_2_plus_v_2 = numpy.power(u, 2) + numpy.power(v, 2)
        u_2_plus_v_2_plus_w_2 = u_2_plus_v_2 + numpy.power(w, 2)

        cos = numpy.cos(rotate)
        sin = numpy.sin(rotate)

        ux_plus_vy_plus_wz = numpy.sum(
            self.__parent_molecule.get_coordinates() * delta, 1
        )

        # Now rotate molecule. In a perform world, I'd have an awesome better
        # numpified version of this, perhaps with tensor or matrix
        # multiplication

        for t in range(len(self.__parent_molecule.get_coordinates())):
            # so t is an atom index
            x_not, y_not, z_not = self.__parent_molecule.get_coordinates()[t]
            ux_plus_vy_plus_wz = u * x_not + v * y_not + w * z_not

            self.__parent_molecule.get_coordinates()[t][0] = (
                a * v_2_plus_w_2 + u * (-b * v - c * w + ux_plus_vy_plus_wz) +
                (
                    -a * v_2_plus_w_2 +
                    u * (b * v + c * w - v * y_not - w * z_not) +
                    v_2_plus_w_2 * x_not
                ) * cos +
                numpy.sqrt(u_2_plus_v_2_plus_w_2) *
                (-c * v + b * w - w * y_not + v * z_not) * sin
            ) # /u_2_plus_v_2_plus_w_2

            self.__parent_molecule.get_coordinates()[t][1] = (
                b * u_2_plus_w_2 + v * (-a * u - c * w + ux_plus_vy_plus_wz) +
                (
                    -b * u_2_plus_w_2 +
                    v * (a * u + c * w - u * x_not - w * z_not) +
                    u_2_plus_w_2 * y_not
                ) * cos +
                numpy.sqrt(u_2_plus_v_2_plus_w_2) *
                (c * u - a * w + w * x_not - u * z_not) * sin
            ) # /u_2_plus_v_2_plus_w_2

            self.__parent_molecule.get_coordinates()[t][2] = (
                c * u_2_plus_v_2 + w * (-a * u - b * v + ux_plus_vy_plus_wz) +
                (
                    -c * u_2_plus_v_2 +
                    w * (a * u + b * v - u * x_not - v * y_not) +
                    u_2_plus_v_2 * z_not
                ) * cos +
                numpy.sqrt(u_2_plus_v_2_plus_w_2) *
                (-b * u + a * v - v * x_not + u * y_not) * sin
            ) # /u_2_plus_v_2_plus_w_2

        self.__parent_molecule.set_coordinates(
            self.__parent_molecule.get_coordinates() *
            (1.0 / u_2_plus_v_2_plus_w_2)
        )

        # here I'm going to just delete the hierarchical info because I'm lazy.
        try:
            # calculated bounding spheres, if any, are no longer valid.
            del self.__parent_molecule.get_hierarchy()['spheres']
        except:
            pass

    def rotate_molecule_around_a_line_between_atoms(self, line_point1_index,
                                                    line_point2_index, rotate):
        """Rotate the molecular model about a line segment. The end points of
        the line segment are atoms of specified indices.

            Args:
                line_point1_index -- An int, the index of the first atom at one
                    end of the line segment.
                line_point2_index -- An int, the index of the second atom at
                    the other end of the line segment.
                rotate -- A float, the angle of rotation, in radians.

            """

        pt1 = self.__parent_molecule.get_coordinates()[line_point1_index]
        pt2 = self.__parent_molecule.get_coordinates()[line_point2_index]
        self.rotate_molecule_around_a_line_between_points(pt1, pt2, rotate)

        try:
            # calculated bounding spheres, if any, are no longer valid.
            del self.__parent_molecule.get_hierarchy()['spheres']
        except:
            pass

    def rotate_molecule_around_pivot_point(self, pivot, thetax,
                                           thetay, thetaz):
        """Rotate the molecular model around a specified atom.

            Args:
                pivot -- A numpy.array, the (x, y, z) coordinate about which
                    the molecular model will be rotated.
                thetax -- A float, the angle to rotate relative to the x axis,
                    in radians.
                thetay -- A float, the angle to rotate relative to the y axis,
                    in radians.
                thetaz -- A float, the angle to rotate relative to the z axis,
                    in radians.

            """

        if pivot.shape == (3,): pivot = numpy.array([pivot])

        # First, move the Molecule so the pivot is at the origin
        self.__parent_molecule.set_coordinates(
            self.__parent_molecule.get_coordinates() - pivot
        )

        # do the rotation
        sinx = numpy.sin(thetax)
        siny = numpy.sin(thetay)
        sinz = numpy.sin(thetaz)
        cosx = numpy.cos(thetax)
        cosy = numpy.cos(thetay)
        cosz = numpy.cos(thetaz)

        rot_matrix = numpy.array([
            [(cosy * cosz),
             (sinx * siny * cosz + cosx * sinz),
             (sinx * sinz - cosx * siny * cosz)
            ],
            [-(cosy * sinz),
             (cosx * cosz - sinx * siny * sinz),
             (cosx * siny * sinz + sinx * cosz)
            ],
            [siny, -(sinx * cosy), (cosx * cosy)]
        ])

        self.__parent_molecule.set_coordinates(
            numpy.dot(rot_matrix, self.__parent_molecule.get_coordinates().T).T
        )

        # now move the pivot point back to it's old location
        self.__parent_molecule.set_coordinates(
            self.__parent_molecule.get_coordinates() + pivot
        )

        try:
            # calculated bounding spheres, if any, are no longer valid.
            del self.__parent_molecule.information.hierarchy['spheres']
        except:
            pass

    def rotate_molecule_around_pivot_atom(self, pivot_index, thetax,
                                          thetay, thetaz):
        """Rotate the molecular model around a specified atom.

            Args:
                pivot_index -- An int, the index of the atom about which the
                    molecular model will be rotated.
                thetax -- A float, the angle to rotate relative to the x axis,
                    in radians.
                thetay -- A float, the angle to rotate relative to the y axis,
                    in radians.
                thetaz -- A float, the angle to rotate relative to the z axis,
                    in radians.

            """

        pivot = self.__parent_molecule.get_coordinates()[pivot_index]
        self.rotate_molecule_around_pivot_point(pivot, thetax, thetay, thetaz)

        try:
            # calculated bounding spheres, if any, are no longer valid.
            del self.__parent_molecule.get_hierarchy()['spheres']
        except:
            pass
