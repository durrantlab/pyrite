# Plan:
# Hard code input for spheres
# Find all points in a sphere
# Select all bones in a sphere
# Keep only every certain number of bones in the sphere (delete others)
 
import mathutils
from bpy import context

# Hard code input for spheres
spheres = [
    # [x, y, z, radius, every_other]
    [0, 0, 0, 25, 10],
    [-30, -70, 375, 25, 100],
]

# Create kd-tree containing current bone coordinates
obj = context.object
mesh = obj.data
# size = len(mesh.vertices) # Mesh is an armature, doesn't have vertices attribute
size = len(bpy.context.visible_bones)   # Tree size is current number of bones
kd = mathutils.kdtree.KDTree(size)
bone_list = []

for i in range(size):
    bone_list.append(bpy.context.visible_bones[0].head)

for i, v in enumerate(bone_list):
    kd.insert(v, i)

kd.balance()

# Find all bones within radius
co_find = (-30, -70, 375)
for (co, index, dist) in kd.find_range(co_find, 100):
    print(" hi ", co, index, dist)

# Keep only every certain number of bones
for x in spheres:
    co_find = (x[0], x[1], x[2])
    radius = x[3]
    every_other = x[4]
    print(kd.find_range(co_find, radius))
