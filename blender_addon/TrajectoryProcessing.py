# Pyrite is a Blender addon for visualization molecular dynamics simulations.
# Copyright (C) 2017  Jacob D. Durrant
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import bpy
import mathutils
import numpy
from . import scoria
from .DurBlend import BackgroundJobParentClass
from .DurBlend import Messages
from . import globals

try: import cStringIO as StringIO
except: from io import StringIO

class ProcessTrajectory(BackgroundJobParentClass):
    """
    A class to load and process a molecular dynamics trajectory.
    """
    
    bl_idname = "process.trajectory"

    def setup(self, context, event):
        """
        Set up variables needed to process the trajectory.

        :param bpy_types.Context context: The context.

        :param bpy.types.Event event: The event.
        """
        
        self.overall_pruning_stride = context.object.overall_pruning_stride
        self.current_step = "START"
        self.current_frame = None
        self.pruning_spheres = []
        self.selection_atoms_to_keep_with_offset = {}
        self.selection_atoms_to_keep = []
        self.frame_index = 0
        self.frame_stride = context.object.frame_stride
        self.pdb_filename = context.object.pdb_filename
        self.frames = self.get_frames()
        self.guide_empties_location = None
        self.protein_obj = bpy.context.scene.objects.active

    def run_step(self, context, event):
        """
        Run a single step of the background job. You need to do it in steps so
        you can periodically return control to the UI, making it seem like a
        background job.

        :param bpy_types.Context context: The context. 

        :param bpy.types.Event event: The event.
        """
        
        if self.current_step == "START":
            # On the first step, set the current frame and post a message.
            self.current_frame = next(self.frames)
            Messages.send_message(
                "LOAD_TRAJ_PROGRESS", 
                "Identifying which atoms to keep..."
            )
            self.current_step = "ID_ATOMS"
        elif self.current_step == "ID_ATOMS":
            # Second step. Loading the first frame.
            self.get_pruned_indecies()
            self.add_empties_at_atom_points()
            Messages.send_message("LOAD_TRAJ_PROGRESS", "Loading frame 1...")
            self.current_step = "POSITION_EMPTIES"
        if self.current_step == "LOAD_FRAME":
            # First half if third step. Loading additional frames...
            try:
                self.current_frame = next(self.frames)
                Messages.send_message(
                    "LOAD_TRAJ_PROGRESS", 
                    "Loading frame " + str(self.frame_index + 1) + "..."
                )
                self.current_step = "POSITION_EMPTIES"
                return None
            except StopIteration:
                # Here you've reached the last frame, so proceed to the next
                # step.
                self.position_empties_at_atom_locs(position_all=True)
                Messages.send_message(
                    "LOAD_TRAJ_PROGRESS", 
                    "Setting up armature and bones..."
                )
                self.current_step = "ARMATURE"
        elif self.current_step == "POSITION_EMPTIES":
            # Second half of third step. Position empties at the location of
            # the atoms.
            self.position_empties_at_atom_locs()
            self.frame_index = self.frame_index + 1
            self.current_step = "LOAD_FRAME"
        elif self.current_step == "ARMATURE":
            # Everything loaded. Now set up armature and bones.
            self.armature_and_bones()
            return {'CANCELLED'}
    
    def job_cancelled(self):
        """
        The user cancels the job.
        """

        # global currently_loading_traj
        bpy.ops.remove.animations('INVOKE_DEFAULT')
        globals.currently_loading_traj = False

    def get_frames(self):
        """
        Get the next frame from the trajectory. Because this is a generator,
        we only load in one frame at a time. That's good for memory.
        """

        current_frame = 0
        pdb_file = open(self.pdb_filename, 'r')
        current_pdb_frame_lines = []

        while True:
            line = pdb_file.readline()
            current_pdb_frame_lines.append(line)
            if line.strip() in ["END", "ENDMDL"]:
                print("Frame", current_frame)

                mol = scoria.Molecule()
                mol.load_pdb_into_using_file_object(
                    StringIO("".join(current_pdb_frame_lines)),
                    False, False, False, False
                )
                yield mol   

                # Start over with next frame
                current_pdb_frame_lines = []
                current_frame += 1

            if line == "":  # EOF
                break

    def add_empties_at_atom_points(self):
        """
        Add enough empty objects for the number of atoms that will be loaded.
        """

        # Make an empty that all the other empties will be parented to. This
        # is just to keep things organized.
        guide_empties = bpy.ops.object.empty_add(type='PLAIN_AXES', radius=1)
        guide_empties = bpy.context.object
        guide_empties.name = "Pyrite_GuideEmpties"
        self.guide_empties_location = numpy.mean(
            self.current_frame.get_coordinates(), 
            axis=0
        )
        guide_empties.location = self.guide_empties_location
        
        # Add enough empties to match the number of bones. 
        for index in self.selection_atoms_to_keep:
            empty = bpy.ops.object.empty_add(type='PLAIN_AXES', radius=1)
            empty = bpy.context.object
            empty.name = "Pyrite_empty" + str(index)
            empty.parent = guide_empties

    def get_pruned_indecies(self):
        """
        Applies pruning spheres to the existing protein. Determines the
        indecies of the atoms to keep.
        """

        # The key is to use the smallest pruning stride possible for a given
        # point.

        # Add the pruning spheres
        # Add oveall pruning sphere
        self.pruning_spheres.append((self.overall_pruning_stride, 0.0, 0.0, 0.0, 1e50))

        # Add individual spheres
        for sphere in [obj for obj in bpy.data.objects if obj.name.startswith("Pyrite_highres_sphere__")]:
            x, y, z = list(sphere.location)
            r = 5.0 * sphere.scale.x  # Because radius is 5 at creation
            stride = sphere.sphere_pruning_stride
            self.pruning_spheres.append((stride, x, y, z, r))

        # Make a kd tree if needed
        coors = self.current_frame.get_coordinates()  # So kdtree calculated on
                                                      # coordinates of first frame only.
        # Create kd-tree containing atom/bone coordinates
        kdtree = mathutils.kdtree.KDTree(len(coors))
        bone_list = []

        # Add coordinates to tree
        for i, c in enumerate(coors):
            kdtree.insert(c, i)
        kdtree.balance()

        # Make sure the pruning spheres are ordered by the stride, from
        # smallest stride to greatest.
        self.pruning_spheres.sort()

        # Go through each sphere and apply a mask, where true means the
        # coordinate is in the given sphere, and false means it isn't.
        indices_in_previous_spheres = set([])
        total_indices_to_keep = set([])
        for sphere in self.pruning_spheres:
            # Get the coordinate indices that are in the sphere
            atom_stride, center_x, center_y, center_z, radius = sphere
            co_find = (center_x, center_y, center_z)
            coors_in_sphere = kdtree.find_range(co_find, radius)
            # coors_in_sphere looks like this:
            # (Vector((-26.695999145507812, -92.79199981689453, 389.52801513671875)), 777, 10.07503890991211)

            # Just keep the indices in this sphere
            indices_in_this_sphere = set([])
            for coor in coors_in_sphere:
                indices_in_this_sphere.add(int(coor[1]))
            
            # Remove from the ones in this sphere ones that were in previous spheres
            indices_to_keep = indices_in_this_sphere - indices_in_previous_spheres

            # Update indices_in_previous_spheres
            indices_in_previous_spheres = indices_in_previous_spheres.union(indices_to_keep)

            # Now only keep every few of the indices_to_keep
            indices_to_keep_sparce = numpy.array(list(indices_to_keep))[::atom_stride]

            # Update all the totals
            total_indices_to_keep = total_indices_to_keep.union(set(indices_to_keep_sparce))

        # Keep track of the atomic indecies to keep.
        total_indices_to_keep = list(total_indices_to_keep)
        total_indices_to_keep.sort()
        total_indices_to_keep = numpy.array(total_indices_to_keep)
        self.selection_atoms_to_keep = total_indices_to_keep
        for i in range(self.frame_stride):
            self.selection_atoms_to_keep_with_offset[i] = total_indices_to_keep[i::self.frame_stride]

    def position_empties_at_atom_locs(self, position_all=False):
        """
        Positions the empties at the locations of the retained atoms.

        :param boolean position_all: Whether or not to position all empties.
                       Use True at first and last frame. Defaults to False.
        """

        try:  
            # So dumb that blender throws an error if it's already in object
            # mode...
            bpy.ops.object.mode_set(mode='OBJECT')
        except:
            pass

        # Sets next frame to add
        bpy.context.scene.frame_set(self.frame_index)

        if self.frame_index == 0 or position_all == True:
            # Set all of them.
            sel = self.selection_atoms_to_keep
        else:
            # Set just the ones appropriate for this frame (using offset).
            sel = self.selection_atoms_to_keep_with_offset[
                self.frame_index % self.frame_stride
            ]
        
        # Get the coordinates of the appropriate atoms.
        coors = self.current_frame.get_coordinates()[sel]

        # Go through and set the appropriate empty object to the associated
        # coordinate.
        for coor_index, coor in enumerate(coors):
            si = sel[coor_index]
            empty = bpy.data.objects["Pyrite_empty" + str(si)]
            empty.location = coor - self.guide_empties_location
            empty.keyframe_insert(data_path='location')

    def armature_and_bones(self):
        """
        Setup the armature and bones.
        """

        # global currently_loading_traj

        # Creating armature
        bpy.ops.object.add(type='ARMATURE', enter_editmode=True)
        object = bpy.context.object
        object.name = 'Pyrite_Armature'
        armature = object.data
        armature.name = 'Frame'

        # Add bones
        for index in self.selection_atoms_to_keep:
            bone_name = 'Pyrite_bone' + str(index)
            bone = armature.edit_bones.new(bone_name)
            bone.head = (0, 0, 0)
            bone.tail = (0, 0, 2)
            bone.envelope_weight = 1.0  # Needed for envelope-based mesh
                                        # vertex weighting.
            bone.envelope_distance = 5.0

        # Now constrain them to the locations of the empty objects (which are
        # animated).
        bpy.ops.object.mode_set(mode='POSE')
        armature = bpy.data.objects["Pyrite_Armature"]

        for index in self.selection_atoms_to_keep:
            bone = armature.pose.bones['Pyrite_bone' + str(index)]
            constraint = bone.constraints.new(type="COPY_LOCATION")
            constraint.target = bpy.data.objects["Pyrite_empty" + str(index)]

        # Now make sure the pose at frame 0 is set as the rest pose (so you
        # can do automatic weights later...)
        bpy.context.scene.frame_set(0)
        bpy.ops.pose.armature_apply()

        # Parent the protein mesh to that armature
        for obj in bpy.data.objects: obj.select = False
        self.protein_obj.select = True
        armature.select = True
        bpy.context.scene.objects.active = armature
        bpy.ops.object.parent_set(type="ARMATURE_AUTO")

        # Zoom in on armature
        try: bpy.ops.object.mode_set(mode='OBJECT')
        except: pass
        
        bpy.context.scene.objects.active = bpy.data.objects['Pyrite_Armature']
        bpy.ops.view3d.view_selected(use_all_regions=False)

        # No longer running
        globals.currently_loading_traj = False
