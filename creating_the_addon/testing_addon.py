try: import cStringIO as StringIO
except: from io import StringIO

import scoria
import numpy
import json
import mathutils
import bpy
from bpy import context
from bpy.props import *

bl_info = {
    "name": "Mineral",
    "author" : "Name <name@example.com>",
    "version" : (1, 0, 0),
    "blender" : (2, 5, 7),
    "location" : "View 3D > ",
    "description" : "",
    "warning" : "",
    "wiki_url" : "",
    "tracker_url" : "",
    "category": "Object",
}

##### Setup scene and object variables #####
def nothing(self, context):
    """
    This is a function that does nothing.
    """
    return

class Properties:
    """
    This class contains functions that easily define the properties controlled
    by a given widget. These properties can be seen as typed variables with
    associated descriptions. For example, suppose you want to have a checkbox
    in your widget to indicate whether or not fainting goats are awesome. That
    would be a boolean, so the property could be defined using the boolProp
    definition in this class.
    """

    def intProp(self, txt, min=-100, max=100, default=33, update=nothing):
        """
        Define an integer property.

        :param str txt: A prompt to let the user know what the property is
                   for. For example, "How many fainting goats are there?"

        :param int min: The minimum value this property can have. Defaults to
                   -100.

        :param int min: The maximum value this property can have. Defaults to
                   100.

        :param int default: The default value of this property. Defaults
                   to 33.

        :param func update: The function to execute when this value is
                    updated. Defaults to the nothing function defined above,
                    which does nothing.

        :returns: a dictionary with the specified values.
        :rtype: :class:`str`  # What is this line?
        """

        return IntProperty(
            name=txt,
            min=min, max=max,
            default=default,
            description="An integer between " + str(min) + " and " + str(max),
            update=update
        )

    # MAKE FLOAT VECTOR PROPERTY
    # def intVectorProp(self, txt, min=(-1000, -1000, -1000), max=(1000, 1000, 1000), default=(0, 0, 0), subtype='NONE', size=3, update=nothing):
    #     """
    #     Define an integer vector property.

    #     :param str txt: A prompt to let the user know what the property is
    #                for. For example, "How many fainting goats are there?"

    #     :param int min: The minimum value this property can have. Defaults to
    #                -100.

    #     :param int min: The maximum value this property can have. Defaults to
    #                100.

    #     :param int default: The default value of this property. Defaults
    #                to 33.

    #     :param func update: The function to execute when this value is
    #                 updated. Defaults to the nothing function defined above,
    #                 which does nothing.

    #     :returns: a dictionary with the specified values.
    #     :rtype: :class:`str`  # What is this line?
    #     """

    #     return IntVectorProperty(
    #         name=txt,
    #         min=min, max=max,
    #         default=default,
    #         subtype=subtype,
    #         size=size,
    #         description="A vector between " + str(min) + " and " + str(max),
    #         update=update
    #     )

    def floatProp(self, txt, min=-100.0, max=100.0, default=33.0, update=nothing):
        """
        Define a float property.

        :param str txt: A prompt to let the user know what the property is
                   for. For example, "How many fainting goats are there,
                   including fractional goats?"

        :param int min: The minimum value this property can have. Defaults to
                   -100.0.

        :param int min: The maximum value this property can have. Defaults to
                   100.0.

        :param int default: The default value of this property. Defaults
                   to 33.0.

        :param func update: The function to execute when this value is
                    updated. Defaults to the nothing function defined above,
                    which does nothing.

        :returns: a dictionary with the specified values.
        :rtype: :class:`str`  # What is this line?
        """

        return FloatProperty(
            name=txt,
            min=min, max=max,
            default=default,
            description="A float between " + str(min) + " and " + str(max),
            update=update
        )

    def boolProp(self, txt, default=True, update=nothing):
        """
        Define a boolean property.

        :param str txt: A prompt to let the user know what the property is
                   for. For example, "Are fainting goats great?"

        :param bool default: The default value of this property. Defaults
                    to True.

        :param func update: The function to execute when this value is
                    updated. Defaults to the nothing function defined above,
                    which does nothing.

        :returns: a dictionary with the specified values.
        :rtype: :class:`str`  # What is this line?
        """

        return BoolProperty(
            name=txt,
            default=default,
            description="True or false",
            update=update
        )

    def strProp(self, txt, default="", subtype='NONE', update=nothing):
        """
        Define a string property.

        :param str txt: A prompt to let the user know what the property is
                   for. For example, "Type in the name of your favorite
                   fainting goat."

        :param str default: The default value of this property. Defaults
                    to "".

        :param func update: The function to execute when this value is
                    updated. Defaults to the nothing function defined above,
                    which does nothing.

        :returns: a dictionary with the specified values.
        :rtype: :class:`str`  # What is this line?
        """

        return StringProperty(
            name=txt,
            default=default,
            description="Text",
            subtype=subtype,
            update=update
        )

    def enumProp(self, txt, items=[("moose", "Moose", ""), ("dog", "Dog", "")], update=nothing):
        """
        Define an enumerated property.

        :param str txt: A prompt to let the user know what the property is
                   for. For example, "Which of these is a king of goat?"

        :param ??? items: A list of tuples. Each tuple represents an option.
                   The first item in the tuple is the option name. The second
                   item is the option name in a more human-readable format.
                   The third item is the value if this option is selected.

        :param func update: The function to execute when this value is
                    updated. Defaults to the nothing function defined above,
                    which does nothing.

        :returns: a dictionary with the specified values.
        :rtype: :class:`str`  # What is this line?
        """

        return EnumProperty(
            name=txt,
            #default = items[0],
            description="Select Option",
            update=update,
            items=items
        )

##### Class for drawing UI elements #####
class UI:
    """
    This class contains functions to make it easier to layout the user
    interface of a blender addon panel.
    """

    row_context = None
    parent = None

    def use_layout_row(self):
        """
        Tells the UI to use a layout row rather than a box row. The following
        widgets (rows) are not grouped in a box, and a group label isn't added
        by default.
        """

        self.row_context = self.parent.layout

    def use_box_row(self, label_txt):
        """
        Tells the UI to use a box row rather than a layout row. The widgets
        (rows) that follow are grouped in a box. A group title is added.

        :param str label_txt: The group title.
        """

        box = self.parent.layout.box()
        box.label(label_txt)
        self.row_context = box

    def new_row(self):
        """
        Start a new row. Whether it is a layout row or a box row depends on
        whether use_box_row() or use_layout_row() was called (above).
        """

        row = self.row_context.row(align=True)
        row.alignment = "EXPAND"
        return row

    def label(self, txt="Label Text"):
        """
        Add a simple text label to the current row.

        :param str txt: The label. Defaults to "Label Text".
        """

        row = self.new_row()
        row.label(text=txt)

    def object_property(self, property_name="location"):
        """
        Add a widget to the addon that controls a given object. How do you
        know if this widget will serve an intger, float, boolean, etc.
        property? That is specified in the associated property data, set using
        the functions of the Properties class.

        :param str property_name: The name of the property. With this name,
                   the code will look up the original Property type and will
                   add the appropriate widget to your addon.
        """

        row = self.new_row()
        row.prop(self.parent.obj, property_name)

    def scene_property(self, property_name="location"):
        """
        Add a widget to the addon that controls some aspect of the entire
        scene. How do you know if this widget will serve an intger, float,
        boolean, etc. property? That is specified in the associated property
        data, set using the functions of the Properties class.

        :param str property_name: The name of the property. With this name,
                   the code will look up the original Property type and will
                   add the appropriate widget to your addon.
        """

        row = self.new_row()
        row.prop(self.parent.scene, property_name)

    def ops_button(self, rel_data_path="object.modifier_add", button_label="Add Modifier!"):
        """
        Add a button to your widget. Use this button when you want to do
        something simple, without passing any additional actions (which are
        like parameters to get a given functionality to do something more
        specific that what is generic).

        So, for example, if you want to select/deselect all the objects in
        your scene, you could use this button with the rel_data_path set to
        "object.select_all".

        :param str rel_data_path: A string specifiying what this button should
                   do.

        :param str button_label: The text of the button.
        """

        # Note that rel_data_path does not include bpy.ops.
        # So instead of bpy.ops.object.modifier_add, just object.modifier_add
        row = self.new_row()
        row.operator(rel_data_path, text=button_label) #, icon='FILESEL')

    def ops_action_button(self, rel_data_path="object.select_all", button_label="Invert Selection!", action="INVERT"):
        """
        Add a button to your widget. Use this button when you want to do
        something beyond generic functionality (i.e., when you need to use an
        "action", which is like a parameter to get a given functionality to do
        something more specific that what is generic).

        So, for example, if you want to invert the current selection, setting
        the rel_data_path set to "object.select_all" won't do. You
        additionally need to use the "INVERT" action.

        :param str rel_data_path: A string specifiying what this button should
                   do.

        :param str button_label: The text of the button.

        :param str action: The action to perform when the button is pressed.
        """

        row = self.new_row()
        row.operator(rel_data_path, text=button_label).action = action

class PanelParentClass(bpy.types.Panel):
    """
    This class is the parent class of any widget-specific panel class you
    might make. Don't change this class, but change a class of your own that
    inherits this one.
    """

    # All panels will have associated objects and an associated scene. Obj is
    # the last object selected. So obj.name, for example, is its name. The
    # value shown in your widget will update automatically (you don't need to
    # explicitly draw it yourself because the draw function is called
    # frequently).
    obj = None
    scene = None

    # All panels will have associated properties and a user interface, so make
    # those here.
    prop_funcs = Properties()
    ui = UI()

    @classmethod
    def start(self):
        """
        This function is called when your panel is created. Every panel must
        have one. In this case, it just calls setup_properties(), which you
        define in your Panel class. I'm keeping these separate in case in the
        future we need to make sure some code is run when initializing any
        panels.

        Note that setup_properties() must be a classmethod. It's a place where
        you define all the properties for your panel (see example below).
        """

        self.setup_properties()

    @classmethod
    def setup_properties(self):
        """
        This function should be overwritten in your child class. It's a place
        where you define all the properties for your panel (see example
        below).
        """

        assert False, "You need to define a setup_properties() definition in your own Panel class!"

    def set_class_variables(self, context):
        self.obj = context.object
        self.scene = bpy.context.scene
        self.ui.parent = self



class Mineral(PanelParentClass):
    """Mineral"""
    bl_label = "Mineral"
    bl_idname = "object.mineral"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'

    @classmethod
    def setup_properties(self):
        """
        Define all the scene and object properties for your panel. Every Panel
        class must have this function!
        """

        # frame_stride description = Every # of frames to keep
        # overall_pruning_stride description = Every # of atoms to keep

        # Set up scene and object properties.
        bpy.types.Object.pdb_filename = self.prop_funcs.strProp("Select PDB file", "sample.pdb", 'FILE_PATH', nothing)
        bpy.types.Object.frame_stride = self.prop_funcs.intProp("Frame stride", 1, 100, 2, nothing)
        
        bpy.types.Object.overall_pruning_stride = self.prop_funcs.intProp("Overall atom stride", 1, 100, 5, nothing)

        bpy.types.Object.sphere_coordinate = bpy.props.FloatVectorProperty(
            name="Center coordinates",
            default=(0.0, 0.0, 0.0)
        )
        bpy.types.Object.sphere_radius = self.prop_funcs.intProp("Radius of sphere", 1, 100, 20, nothing)
        bpy.types.Object.sphere_pruning_stride = self.prop_funcs.intProp("Pruning stride", 1, 100, 5, nothing)

        # bpy.types.Object.center_coord = self.prop_funcs.intVectorProp("Center coordinates", (-1000, -1000, -1000), (1000, 1000, 1000), (0, 0, 0), 'NONE', 3, nothing)
        # Need a float vector property to input atom coordinates
        # How do we make it possible to enter more coordinates if desired?

    def draw(self, context):
        """
        Every panel class must have a draw function. It gets called over and
        over again, constantly refreshing the Panel's appearance (and updating
        the displayed values in the process).

        THE DRAW FUNCTION MUST ALWAYS START WITH:
        self.set_class_variables(context)

        :param ??? context: The context of the currently selected object.
        """

        self.set_class_variables(context)

        self.ui.use_layout_row()
        self.ui.label("Protein (Object Name: " + self.obj.name + ")")

        self.ui.use_box_row("Load a Protein")
        self.ui.object_property(property_name="pdb_filename")
        self.ui.object_property(property_name="frame_stride")
        self.ui.new_row()

        self.ui.use_box_row("Pruning")
        self.ui.object_property(property_name="overall_pruning_stride")
        self.ui.new_row()

        self.ui.use_layout_row()
        self.layout.operator("protein.display")

        self.ui.new_row()
        self.ui.new_row()

        self.ui.use_box_row("Add Pruning Sphere")
        self.ui.object_property(property_name="sphere_coordinate")
        self.ui.object_property(property_name="sphere_radius")
        self.ui.object_property(property_name="sphere_pruning_stride")

        self.layout.operator("add.sphere")

        self.ui.new_row()
        self.layout.operator("finalize.spheres")

def load_pdb_trajectory(pdb_filename, frame_stride):
    """
    Loads molecule trajectory from a given PDB file into a numpy array.
    Keeps every certain number of frames in array based on user-inputted value.

    Args:
    pdb_filename (string): The name of a PDB file (including '.pdb').
    frame_stride (integer): The stride for frames to keep.

    Returns:

    """
    obj = context.object

    print("Loading PDB trajectory: " + pdb_filename)
    pdb_filename = pdb_filename
    frame_stride = frame_stride

    # Load the trajectory
    trajectory = scoria.Molecule()
    trajectory.load_pdb_trajectory_into(pdb_filename, bonds_by_distance = False, serial_reindex = False, resseq_reindex = False)

    # Delete every frame_stride frames.
    print("Keeping only every " + str(frame_stride) + " frames...")
    frame_indices = numpy.array(range(trajectory.get_trajectory_frame_count()))
    frame_indices_to_keep = frame_indices[::frame_stride]
    frame_indices_to_delete = numpy.setdiff1d(frame_indices, frame_indices_to_keep)
    for idx in frame_indices_to_delete[::-1]:
        trajectory.delete_trajectory_frame(idx)
    return trajectory



def add_overall_pruning_stride(pruning_spheres, atom_stride):
    """
    Prunes the protein based on the user-inputted atom stride. Deletes every
    certain number of atoms in the coordinate array based on user-inputted value.

    Args:
    atom_stride (integer): The stride for atoms to keep in the entire protein.

    Returns:

    """
    print("Keeping only every " + str(atom_stride) + " atoms...")
    pruning_spheres.append((atom_stride, 0.0, 0.0, 0.0, 1e50))
    return pruning_spheres

def add_pruning_sphere(pruning_spheres, center_x, center_y, center_z, radius, atom_stride):
    """
    Prunes atoms in the protein within a sphere of a certain radius around set coordinates.
    Atom stride, coordinates, and radius are set by user input.

    Args:
    center_x (integer): x coordinate for center of sphere.
    center_y (integer): y coordinate for center of sphere.
    center_z (integer): z coordinate for center of sphere.
    radius (integer): Radius of the sphere.
    atom_stride (integer): Stride for atoms to keep within the sphere.

    Returns:

    """
    pruning_spheres.append((atom_stride, center_x, center_y, center_z, radius))
    return pruning_spheres

def apply_prune(trajectory, kdtree, pruning_spheres):
    """
    Applies pruning spheres to the existing protein.

    Args:

    Returns:

    """
    # The key is to use the smallest pruning stride possible for a given
    # point.

    # Make a kd tree if needed
    coors = trajectory.get_coordinates(frame=0)  # So kdtree calculated on
                                                 # coordinates of first frame only.
    if kdtree is None:
        # Create kd-tree containing atom/bone coordinates
        kdtree = mathutils.kdtree.KDTree(len(coors))
        bone_list = []

        # Add coordinates to tree
        for i, c in enumerate(coors):
            kdtree.insert(c, i)

        kdtree.balance()

    # Make sure the pruning spheres are ordered by the stride, from
    # smallest to greatest.
    pruning_spheres.sort()

    # Go through each sphere and apply a mask, where true means the
    # coordinate is in the given sphere, and false means it isn't.
    masks = []
    for sphere in pruning_spheres:
        # Get the coordinate indices that are in the sphere
        atom_stride, center_x, center_y, center_z, radius = sphere
        co_find = (center_x, center_y, center_z)
        coors_in_sphere = kdtree.find_range(co_find, radius)
        coors_in_sphere = numpy.array(coors_in_sphere)

        # doesn't work.
        indices_in = set([])
        for coor in coors_in_sphere:
            indices_in.add(int(coor[1]))

        # Make the mask, with those set to true that are within the
        # sphere.
        mask = numpy.zeros(trajectory.get_total_number_of_atoms()).astype(bool)
        indices_in = list(indices_in)

        if len(indices_in) > 0:
            mask[indices_in] = True

        # Save that mask
        masks.append(mask)

    # Now go through each of the points and decide whether or not to keep
    # it.
    indices_to_keep = []
    for coor_index, coor in enumerate(coors):
        # Find the sphere that this point is in with the lowest stride.
        # Use that stride.
        for sphere_index, sphere in enumerate(pruning_spheres):
            if masks[sphere_index][coor_index] == True:
                # The point is in one of the spheres.
                atom_stride = pruning_spheres[sphere_index][0]
                if coor_index % atom_stride == 0:
                    # It matches the stride, so keep it.
                    indices_to_keep.append(coor_index)
                break  # No need to keep looking through the spheres for
                        # this point. You've got your match.

    # Actually prune the molecule.
    trajectory = trajectory.get_molecule_from_selection(indices_to_keep)
    return trajectory

def make_bones_from_molecules(trajectory, frame_stride):
    """
    """
    try:  # So dumb that blender throws an error if it's already in object mode...
        bpy.ops.object.mode_set(mode='OBJECT')
    except:
        pass

    # Add enough empties to match the number of bones. 
    for index in range(trajectory.get_total_number_of_atoms()):
        empty = bpy.ops.object.empty_add(type='PLAIN_AXES', radius=1)
        empty = bpy.context.object
        empty.name = "empty" + str(index)

    # Now go through the frames and position those empties
    for frame_index in range(0, trajectory.get_trajectory_frame_count(), frame_stride):
        bpy.context.scene.frame_set(frame_index)  # Sets next frame to add

        for coor_index, coor in enumerate(trajectory.get_coordinates(frame=frame_index)):
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
    for index in range(trajectory.get_total_number_of_atoms()):
        bone_name = 'bone' + str(index)
        bone = armature.edit_bones.new(bone_name)
        bone.head = (0, 0, 0)
        bone.tail = (0, 0, 2)
        #bone.envelope_weight = 1.0  # Needed for envelope-based mesh vertex weighting.
        #bone.envelope_distance = 2.0

    # Now constrain them.
    bpy.ops.object.mode_set(mode='POSE')
    armature = bpy.data.objects["Armature"]

    for index in range(trajectory.get_total_number_of_atoms()):
        bone = armature.pose.bones['bone' + str(index)]
        constraint = bone.constraints.new(type="COPY_LOCATION")
        constraint.target = bpy.data.objects["empty" + str(index)]

    # Now make sure the pose at frame 0 is set as the rest pose (so you
    # can do automatic weights later...)
    bpy.context.scene.frame_set(0)
    bpy.ops.pose.armature_apply()
    
def menu_func(self, context):
    self.layout.operator(Mineral.bl_idname)

class OBJECT_OT_DisplayButton(bpy.types.Operator):
    """
    Button for displaying basic pruned protein.
    """
    bl_idname = "protein.display"
    bl_label = "Display Protein"

    def __init__(self):
        self.trajectory = None
        self.kdtree = None
        self.overall_pruning_stride = 1
        self.pruning_spheres = []
        self.frame_stride = None

    def execute(self, context):
        """
        What should be run when the display button is pressed.
        """

        obj = context.object
        print(obj)
        self.frame_stride = obj['frame_stride']
        self.overall_pruning_stride = obj['overall_pruning_stride']

        self.trajectory = load_pdb_trajectory(obj['pdb_filename'], self.frame_stride)
        self.pruning_spheres = add_overall_pruning_stride(self.pruning_spheres, self.overall_pruning_stride)
        self.pruning_spheres = add_pruning_sphere(self.pruning_spheres, 0, 0, 0, 30, 20)    # Not working right now
        self.trajectory = apply_prune(self.trajectory, self.kdtree, self.pruning_spheres)
        make_bones_from_molecules(self.trajectory, self.frame_stride)

        try:
            bpy.ops.object.mode_set(mode='OBJECT')
        except:
            pass
        bpy.context.scene.objects.active = bpy.data.objects['Armature']
        bpy.ops.view3d.view_selected(use_all_regions=False)
        return{'FINISHED'}

class OBJECT_OT_AddSphereButton(bpy.types.Operator):
    # """
    # Button for adding a positioning sphere.
    # """
    bl_idname = "add.sphere"
    bl_label = "Add a Sphere"

    def execute(self, context):
        """
        Adds a sphere to the scene.
        """

        obj = context.object
        print(obj)

        try:  # So dumb that blender throws an error if it's already in object mode...
            bpy.ops.object.mode_set(mode='OBJECT')
        except:
            pass

        # bpy.ops.mesh.primitive_uv_sphere_add()  # name "pruning_sphere"
        # add sphere at a default location unless optional user input provided
        bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=16, size=8.0, view_align=False, enter_editmode=False, location=(-30, -85, 397), rotation=(0.0, 0.0, 0.0))
        # be able to change size of sphere
        # make sphere more transparent?
        # menu select and edit a sphere already added
        # command to grab object and position (or select object and press G): bpy.ops.transform.translate()

        return{'FINISHED'}

class OBJECT_OT_FinalizeButton(bpy.types.Operator):
    # """
    # Finalize button for removing all positioning spheres.
    # """
    bl_idname = "finalize.spheres"
    bl_label = "Finalize Pruning Spheres"

    def execute(self, context):
        """
        Removes all positioning spheres.
        """

        # delete all objects beginning with "pruning_sphere"
        return{'FINISHED'}

# store keymaps here to access after registration
addon_keymaps = []

##### Registration functions #####
def register():
    """
    Registers this addon.
    """
    Mineral.start()
    bpy.utils.register_class(Mineral)
    bpy.types.VIEW3D_MT_object.append(menu_func)
    bpy.utils.register_class(OBJECT_OT_DisplayButton)
    bpy.utils.register_class(OBJECT_OT_AddSphereButton)
    bpy.utils.register_class(OBJECT_OT_FinalizeButton)

    # # handle the keymap
    # wm = bpy.context.window_manager
    # km = wm.keyconfigs.addon.keymaps.new(name='Object Mode', space_type='EMPTY')
    # kmi = km.keymap_items.new(Mineral.bl_idname, 'SPACE', 'PRESS', ctrl=True, shift=True)
    # # kmi.properties.total = 4
    # addon_keymaps.append(km)

def unregister():
    """
    Good practice to make it possible to unregister addons.
    """

    bpy.utils.unregister_class(__name__)

    bpy.utils.unregister_class(Mineral)
    bpy.types.VIEW3D_MT_object.remove(menu_func)
    bpy.utils.unregister_class(OBJECT_OT_DisplayButton)
    bpy.utils.unregister_class(OBJECT_OT_AddSphereButton)
    bpy.utils.unregister_class(OBJECT_OT_FinalizeButton)

    # # handle the keymap
    # wm = bpy.context.window_manager
    # for km in addon_keymaps:
    #     wm.keyconfigs.addon.keymaps.remove(km)
    # # clear the list
    # del addon_keymaps[:]

if __name__ == "__main__":
    """
    Start the addon!
    """

    register()
