import json

# Loading JSON data into a list in Blender
with open('/Users/niveditarajendiran/Documents/University of Pittsburgh/Senior Year/Durrant Lab/Lab Files/labcode/part2_sparcify_the_trajectory/converting_to_json/frame_coordinate_list.json') as json_data:
    coord_list = json.load(json_data)

# print(coord_list) # Prints entire list of coordinate data
# print(coord_list[0][0]) # [-30.852, -81.458, 365.055] (Coordinates of frame 0, atom 0)

# Creating armature and adding bones (IN PROGRESS)
bpy.ops.object.add(type='ARMATURE', enter_editmode=True)
object = bpy.context.object
object.name = 'Armature'
armature = object.data
armature.name = 'Frame1'

for x in range(len(coord_list[0])):
    bone = armature.edit_bones.new('bonename' + str(x))
    bone.tail = coord_list[0][x]
    bone.head = coord_list[0][x]
    
# Blender keeps crashing after adding all bones
# Try breaking up into multiple armatures

bpy.ops.object.mode_set(mode='OBJECT')


# vec = mathutils.Vector((0.0, 0.0, 1.0))
# eul = mathutils.Euler((0.0, math.radians(45.0), 0.0), 'XYZ')
# bone.tail = mathutils.Vector([0,0,1])
# x = bpy.data.objects["Cube"]
# x.location
# x.location = ([-30.852, -81.458, 365.055])

# x = bpy.data.armatures["Frame1"].bones['bonename2']
# armature = bpy.data.armatures["Frame1"]
# x.head = [-30.852, -81.458, 365.055]