import json
import numpy as np
import mathutils
from scipy import spatial
from bpy import context

# Hard code input for spheres, [x, y, z, radius, every_other]
spheres = [
    [0, 0, 0, 25, 10],
    [-30, -70, 375, 100, 2],
    [-35, -80, 370, 100, 2],
    [-35, -80, 370, 10, 2],
]

# Load the trajectory
with open('/Users/niveditarajendiran/Documents/University of Pittsburgh/Senior Year/Durrant Lab/Lab Files/labcode/part2_sparcify_the_trajectory/get_first_frame/frame_coordinate_list.json') as json_data:
    coordinates = json.load(json_data)

# Create kd-tree containing atom/bone coordinates
# size = len(coordinates)
# kd = mathutils.kdtree.KDTree(size)
# bone_list = []

# Add coordinates to tree
# for i, c in enumerate(coordinates):
#     kd.insert(c, i)

# kd.balance()

kd = scipy.spatial.cKDTree(coordinates)

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