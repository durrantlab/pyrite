try: import cStringIO as StringIO
except: from io import StringIO

import scoria
import numpy
import json
import mathutils
from bpy import context
import bpy

#import _trkr  # Please leave this intact for now.
#_trkr.trkr()

class Mineral:
    def __init__(self):
        print("Creating mineral object!")
        self.trajectory = None
        self.kdtree = None
        self.overall_pruning_stride = 1
        self.pruning_spheres = []
        self.frame_stride = None

    def add_overall_pruning_stride(self, atom_stride):
        print("Keeping only every " + str(atom_stride) + " atoms...")
        self.pruning_spheres.append((atom_stride, 0.0, 0.0, 0.0, 1e50))
        
    def add_pruning_sphere(self, center_x, center_y, center_z, radius, atom_stride):
        self.pruning_spheres.append((atom_stride, center_x, center_y, center_z, radius))

    def apply_prune(self):
        # The key is to use the smallest pruning stride possible for a given
        # point.

        # Make a kd tree if needed
        coors = self.trajectory.get_coordinates(frame=0)  # So kdtree
                                                          # calculated on
                                                          # coordinates of
                                                          # first frame only.
        if self.kdtree is None:
            # Create kd-tree containing atom/bone coordinates
            self.kdtree = mathutils.kdtree.KDTree(len(coors))
            bone_list = []

            # Add coordinates to tree
            for i, c in enumerate(coors):
                self.kdtree.insert(c, i)

            self.kdtree.balance()
        
        # Make sure the pruning spheres are ordered by the stride, from
        # smallest to greatest.
        self.pruning_spheres.sort()

        # Go through each sphere and apply a mask, where true means the
        # coordinate is in the given sphere, and false means it isn't.
        masks = []
        for sphere in self.pruning_spheres:
            # Get the coordinate indecies that are in the sphere
            atom_stride, center_x, center_y, center_z, radius = sphere
            co_find = (center_x, center_y, center_z)
            coors_in_sphere = self.kdtree.find_range(co_find, radius)
            coors_in_sphere = numpy.array(coors_in_sphere)

            # doesn't work.
            indecies_in = set([])
            for coor in coors_in_sphere:
                indecies_in.add(int(coor[1]))

            # Make the mask, with those set to true that are within the
            # sphere.
            mask = numpy.zeros(self.trajectory.get_total_number_of_atoms()).astype(bool)
            indecies_in = list(indecies_in)

            if len(indecies_in) > 0:
                mask[indecies_in] = True

            # Save that mask
            masks.append(mask)
        
        # Now go through each of the points and decide whether or not to keep
        # it.
        indecies_to_keep = []
        for coor_index, coor in enumerate(coors):
            # Find the sphere that this point is in with the lowest stride.
            # Use that stride.
            for sphere_index, sphere in enumerate(self.pruning_spheres):
                if masks[sphere_index][coor_index] == True:
                    # The point is in one of the spheres.
                    atom_stride = self.pruning_spheres[sphere_index][0]
                    if coor_index % atom_stride == 0:
                        # It matches the stride, so keep it.
                        indecies_to_keep.append(coor_index)
                    break  # No need to keep looking through the spheres for
                           # this point. You've got your match.
        
        # Actually prune the molecule.
        self.trajectory = self.trajectory.get_molecule_from_selection(indecies_to_keep)

    def load_pdb_trajectory(self, pdb_filename, frame_stride = 1): #, atom_stride):
        print("Loading PDB trajectory: " + pdb_filename)
        self.pdb_filename = pdb_filename
        self.frame_stride = frame_stride

        # Load the trajectory
        self.trajectory = scoria.Molecule()
        self.trajectory.load_pdb_trajectory_into(self.pdb_filename, bonds_by_distance = False, serial_reindex = False, resseq_reindex = False)

        # Delete every frame_stride frames.
        print("Keeping only every " + str(frame_stride) + " frames...")
        frame_indecies = numpy.array(range(self.trajectory.get_trajectory_frame_count()))
        frame_indecies_to_keep = frame_indecies[::frame_stride]
        frame_indecies_to_delete = numpy.setdiff1d(frame_indecies, frame_indecies_to_keep)
        for idx in frame_indecies_to_delete[::-1]:
            self.trajectory.delete_trajectory_frame(idx)

    def make_bones_from_molecues(self):
        try:  # So dumb that blender throws an error if it's already in object mode...
            bpy.ops.object.mode_set(mode='OBJECT')
        except:
            pass

        # Add enough empties to match the number of bones. 
        for index in range(self.trajectory.get_total_number_of_atoms()):
            empty = bpy.ops.object.empty_add(type='PLAIN_AXES', radius=1)
            empty = bpy.context.object
            empty.name = "empty" + str(index)

        # Now go through the frames and position those empties
        for frame_index in range(0, self.trajectory.get_trajectory_frame_count(), self.frame_stride):
            bpy.context.scene.frame_set(frame_index)  # Sets next frame to add

            for coor_index, coor in enumerate(self.trajectory.get_coordinates(frame=frame_index)):
                empty = bpy.data.objects["empty" + str(coor_index)]
                empty.location = coor
                empty.keyframe_insert(data_path='location')

        # Creating armature
        bpy.ops.object.add(type='ARMATURE', enter_editmode=True)
        object = bpy.context.object
        object.name = 'Armature'
        armature = object.data
        armature.name = 'Frame'

        # Add bones
        for index in range(self.trajectory.get_total_number_of_atoms()):
            bone_name = 'bone' + str(index)
            bone = armature.edit_bones.new(bone_name)
            bone.head = (0,0,0)
            bone.tail = (0,0,2)
            #bone.envelope_weight = 1.0  # Needed for envelope-based mesh vertex weighting.
            #bone.envelope_distance = 2.0
        
        # Now constrain them.
        bpy.ops.object.mode_set(mode='POSE')
        armature = bpy.data.objects["Armature"]

        for index in range(self.trajectory.get_total_number_of_atoms()):
            bone = armature.pose.bones['bone' + str(index)]
            constraint = bone.constraints.new(type="COPY_LOCATION")
            constraint.target = bpy.data.objects["empty" + str(index)]
        
        # Now make sure the pose at frame 0 is set as the rest pose (so you
        # can do automatic weights later...)
        bpy.context.scene.frame_set(0)
        bpy.ops.pose.armature_apply()

    
mineral = Mineral()
mineral.load_pdb_trajectory("M2_traj.pdb", 2)
mineral.add_overall_pruning_stride(15)
mineral.add_pruning_sphere(-30, -85, 397, 8, 1)
mineral.apply_prune()
mineral.make_bones_from_molecues()
