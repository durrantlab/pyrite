import scoria
import numpy as np
import mathutils
from bpy import context

# For testing, in blender remember to:
# import sys
# sys.path.append("/Users/jdurrant/Documents/Work/durrant_git/nivedita/part2_sparcify_the_trajectory/")
# import os
# os.chdir("/Users/jdurrant/Documents/Work/durrant_git/nivedita/part2_sparcify_the_trajectory/")
# from imp import reload
# Then:
# import protein_pruner
# reload("protein_pruner")

# Hard code input for spheres, [x, y, z, radius, every_other]
spheres = [
    [0, 0, 0, 25, 10],
    [-30, -70, 375, 100, 2],
    [-35, -80, 370, 100, 2],
    [-35, -80, 370, 10, 5],
]

# Load the trajectory
traj = scoria.Molecule()
traj.load_pdb_trajectory_into("M2_traj.pdb")
# traj.load_pdb_trajectory_into("/Users/niveditarajendiran/Documents/University of Pittsburgh/Senior Year/Durrant Lab/Lab Files/labcode/part2_sparcify_the_trajectory/M2_traj.pdb", bonds_by_distance = False, serial_reindex = False, resseq_reindex = False)


coordinates = traj.get_coordinates()    # Stored as numpy array

# Create kd-tree containing atom/bone coordinates
size = len(coordinates)
kd = mathutils.kdtree.KDTree(size)
bone_list = []

# Add coordinates to tree
for i, c in enumerate(coordinates):
    kd.insert(c, i)

kd.balance()

new_spheres = []    # Will contain lists of coordinates for all pruned spheres

# Find all bones within radius, keep only every certain number of bones
for x in spheres:
    co_find = (x[0], x[1], x[2])
    radius = x[3]
    every_other = x[4]
    bones_in_sphere = kd.find_range(co_find, radius)
    bones_in_sphere = np.array(bones_in_sphere)
    bones_to_keep = bones_in_sphere[::every_other]
    new_spheres.append(bones_to_keep)

print(bones_in_sphere)
print("---------")
print(bones_to_keep)
