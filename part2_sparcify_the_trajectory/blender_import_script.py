import json

# Loading JSON data into a list in Blender
with open('/Users/niveditarajendiran/Documents/University of Pittsburgh/Senior Year/Durrant Lab/Lab Files/labcode/part2_sparcify_the_trajectory/converting_to_json/frame_coordinate_list.json') as json_data:
    coord_list = json.load(json_data)

# Creating armature and adding bones
bpy.ops.object.add(type='ARMATURE', enter_editmode=True)
object = bpy.context.object
object.name = 'Armatures'
armature = object.data
armature.name = 'Frame1'

# Adds bones for frame 0
for x in range(len(coord_list[0])):
    bone = armature.edit_bones.new('bonename' + str(x))
    bone.head = coord_list[0][x]
    tail = coord_list[0][x]
    tail[2] = tail[2] - 1   # Reduces Z coordinate by 1 for tail to create bone length
    bone.tail = tail    # Sets tail as modified atom coordinate

bpy.ops.object.mode_set(mode='OBJECT')