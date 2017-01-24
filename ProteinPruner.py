# Plan:
# Hard code input for spheres
# Find all points in a sphere
# Select all bones in a sphere
# Keep only every certain number of bones in the sphere (delete others)
 
import mathutils
from bpy import context

# Hard code input for spheres, [x, y, z, radius, every_other]
spheres = [
    [0, 0, 0, 25, 10],
    [-30, -70, 375, 100, 2],
    [-35, -80, 370, 100, 2],
    [-35, -80, 370, 10, 2],
]

# Create kd-tree containing current bone coordinates
size = len(bpy.context.visible_bones)
kd = mathutils.kdtree.KDTree(size)
bone_list = []

# Compiles list of current bone coordinates
for i in range(size):
    bone_list.append(bpy.context.visible_bones[0].head)

# Adds coordinates to tree
for i, v in enumerate(bone_list):
    kd.insert(v, i)

kd.balance()

# Find all bones within radius, keep only every certain number of bones
for x in spheres:
    co_find = (x[0], x[1], x[2])
    radius = x[3]
    every_other = x[4]

    # print(kd.find_range(co_find, radius))
    bones_in_sphere = kd.find_range(co_find, radius)
    bones_to_prune = []

    # Compile list of bones to prune
    for i, b in enumerate(bones_in_sphere):
        if (i % every_other) != 0:
            bones_to_prune.append(b)
    
    # Select and remove bones in bones_to_prune
    for i, b in enumerate(bones_to_prune):
        
    # bpy.context.visible_bones[0].name
