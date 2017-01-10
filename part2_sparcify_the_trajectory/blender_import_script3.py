import json
import bpy

# Note: Should change frame numbers for pdb files to be same as keyframes

# Loading JSON data into a list in Blender
with open('/Users/niveditarajendiran/Documents/University of Pittsburgh/Senior Year/Durrant Lab/Lab Files/labcode/part2_sparcify_the_trajectory/converting_to_json/frame_coordinate_list.json') as json_data:
    coord_list = json.load(json_data)

# Creating armature
bpy.ops.object.add(type='ARMATURE', enter_editmode=True)
object = bpy.context.object
object.name = 'Armature'
armature = object.data
armature.name = 'Frame'

# Adding bones for frame 0 (original position)
for atom in range(len(coord_list[0])):
    bone = armature.edit_bones.new('atom' + str(atom))
    bone.head = coord_list[0][atom]
    tail = coord_list[0][atom]
    tail[2] = tail[2] - 1   # Reduces Z coordinate by 1 for tail to create bone length
    bone.tail = tail    # Sets tail as modified atom coordinate
    bone.envelope_weight = 1.0  # Needed for envelope-based mesh vertex weighting.
    bone.envelope_distance = 10.0

bpy.ops.object.mode_set(mode='OBJECT')

# Posing bones for all frames
arm = bpy.context.active_object # Armature
bpy.ops.object.mode_set(mode='POSE')

for atom in range(0, len(coord_list[0])):
    bone = arm.pose.bones['atom' + str(atom)]
    keyframe_num = 0    # Resets frame number
    for frame in range(0, len(coord_list)):
        bpy.context.scene.frame_set(keyframe_num)  # Sets next frame to add
        bone.location = coord_list[frame][atom]
        bone.keyframe_insert(data_path="location")

        keyframe_num += 5 # Set as frame_stride (from sparcify.py)
