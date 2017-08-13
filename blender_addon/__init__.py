# Mineral is a Blender addon for visualization molecular dynamics simulations.
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
from bpy import context
from mathutils import Vector
import os

from .DurBlend import Properties
from .DurBlend import UI
from .DurBlend import PanelParentClass
from .DurBlend import ButtonParentClass
from .DurBlend import Messages
from .TrajectoryProcessing import ProcessTrajectory

bl_info = {
    "name": "Mineral",
    "author" : "Name <name@example.com>",
    "version" : (1, 0, 0),
    "blender" : (2, 5, 7),
    "location" : "View 3D > ",
    "description" : "Mineral plugin",
    "warning" : "",
    "wiki_url" : "",
    "tracker_url" : "",
    "category": "Object",
}

###### Below specific to this plugin ######
currently_loading_traj = False


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

        # Set up scene and object properties.
        bpy.types.Object.pdb_filename = self.prop_funcs.strProp("PDB file", "sample.pdb", 'FILE_PATH')
        bpy.types.Object.frame_stride = self.prop_funcs.intProp("Keep every n frames", 1, 100, 2)
        bpy.types.Object.overall_pruning_stride = self.prop_funcs.intProp("Keep every n atoms", 1, 100, 5)
        bpy.types.Object.sphere_pruning_stride = self.prop_funcs.intProp("Keep every n atoms", 1, 100, 2)

    def draw(self, context):
        """
        Every panel class must have a draw function. It gets called over and
        over again, constantly refreshing the Panel's appearance (and updating
        the displayed values in the process).

        THE DRAW FUNCTION MUST ALWAYS START WITH:
        self.set_class_variables(context)

        :param bpy_types.Context context: The context.
        """

        global currently_loading_traj
        
        self.set_class_variables(context)

        # Pick the object. Use active object if self.obj is None.
        obj_to_use = self.obj
        if obj_to_use is None:
            obj_to_use = bpy.context.scene.objects.active
        
        # Consider the possibility that a trajectory is current being
        # loaded... If so, panel should only indicate progress.
        if currently_loading_traj:
            self.draw_loading_trajectory()
            return

        # Consider possibility that nothing is selected/active. You need to
        # provide instructions re. how to proceed.
        mesh_not_selected = False
        if obj_to_use is None:
            # Nothing active
            mesh_not_selected = True

        if mesh_not_selected == False:
            # What is more than one thing is selected? That doesn't work
            # either...
            currently_selected = [
                obj for obj in bpy.data.objects 
                if obj.select == True
            ]
            if len(currently_selected) != 1:
                mesh_not_selected = True

        if mesh_not_selected == True:
            # So no mesh is selected. Provide instructions re. how to proceed.
            self.ui.use_box_row("Instructions")
            
            if bpy.context.scene.objects.active == None:
                # No active object, so they need to select one.
                self.ui.label("Select an object for additional options.")
            else:
                # There is an active object, but still need to select one.
                # Maybe multiple are selected. Since active object available,
                # can also give options for starting over.
                self.draw_no_object_selected_start_over()
            
            # Show them what to cite.
            self.draw_citation()

            return

        # What object name to use in the panel title?
        obj_to_use_name = (
            obj_to_use.name   
            if not obj_to_use.name.startswith("Mineral_highres_sphere__")
            else obj_to_use.name.split("__")[1]
        )

        if obj_to_use.name.startswith("Mineral_highres_sphere__"):
            # It's one of the selection spheres... Provide info/options about
            # that to the user.
            self.draw_high_detail_sphere_panel()
        else:
            # It's not one of the selection spheres. Must be a protein mesh.
            self.draw_protein_mesh_header(obj_to_use_name)

            # Check if the location of the object is ok.
            loc = [round(v, 1) for v in list(obj_to_use.location)]
            rot = [round(v, 1) for v in list(obj_to_use.rotation_euler)]
            scale = [round(v, 1) for v in list(obj_to_use.scale)]
            if loc != [0.0, 0.0, 0.0] or rot != [0.0, 0.0, 0.0] or scale != [1.0, 1.0, 1.0]:
                # The selected mesh must not have location, rotation, and scale at rest.
                self.draw_object_coordinates_not_set_error(loc, rot, scale)
            else:
                # The location, rotation, and scaling are ok, so show normal UI
                self.draw_main_protein_mesh_panel()

    def draw_loading_trajectory(self):
        """
        Panel contents when loading a trajectory.
        """

        self.ui.use_box_row("Loading Trajectory")
        Messages.display_message("LOAD_TRAJ_PROGRESS", self)
        self.ui.label("Press Esc to stop loading...")

    def draw_no_object_selected_start_over(self):
        """
        Panel contents when no object is selected. Gives instructions, and
        gives options of starting over.
        """

        self.ui.label("Select protein object in 3D viewer.")

        # Does a previous run exist? If so, provide the option to
        # start over.
        previous_run_exists = False
        for obj in bpy.data.objects:
            if obj.name.startswith("Mineral_"):
                previous_run_exists = True
                break

        if previous_run_exists:
            self.ui.use_box_row("Previous Runs")
            self.ui.ops_button(
                rel_data_path="remove.animations", 
                button_label="Remove Animations"
            )
            self.ui.ops_button(
                rel_data_path="start.over", 
                button_label="Start Over"
            )

    def draw_citation(self):
        """
        Show the citation in the panel.
        """

        # Make sure they know what to cite!
        self.ui.use_box_row("Citation")
        self.ui.label("If you use Mineral, please cite:")
        self.ui.label("{FULL CITATION HERE}")

    def draw_high_detail_sphere_panel(self):
        """
        The panel to display when a high-detail sphere is selected.
        """

        # It's one of the selection spheres... Provide info/options about
        # that to the user.
        self.ui.use_layout_row()
        self.ui.label("High-Detail Region")
        self.ui.use_box_row("Properties")
        self.ui.label("Move/scale the sphere to encompass the region.")
        self.ui.object_property(property_name="sphere_pruning_stride")
        Messages.display_message("SPHERE_STRIDE_TOO_HIGH", self)
        self.ui.use_box_row("Finalize")
        self.ui.ops_button(
            rel_data_path="backto.protein", 
            button_label="Back to Protein Mesh"
        )
        self.ui.ops_button(
            rel_data_path="delete.region", 
            button_label="Delete Region"
        )

    def draw_protein_mesh_header(self, obj_to_use_name):
        """
        Draw the header of the main protein-mesh panel.

        :param str obj_to_use_name: The object name to display.
        """
        
        # Show the name
        self.ui.use_layout_row()
        self.ui.label("Protein Mesh (Object Name: " + obj_to_use_name  + ")")

        # Provide button to return to the main menu.
        self.ui.ops_button(
            rel_data_path="main.menu", 
            button_label="Return to Main Menu"
        )

    def draw_object_coordinates_not_set_error(self, loc, rot, scale):
        """
        Draw the error panel if the loc, rot, and scale aren't "at rest."

        :param [float, float, float] loc: A list of floats representing the
                                     object location.

        :param [float, float, float] rot: A list of floats representing the
                                     object rotation.

        :param [float, float, float] scale: A list of floats representing the
                                     object scaling.
        """

        # The selected mesh must not have location, rotation, and scale at rest.
        self.ui.use_box_row("Trajectory and Protein-Mesh Positions Must Match!")
        if loc != [0.0, 0.0, 0.0]:
            self.ui.label("Mesh location " + str(loc) + " is not [0.0, 0.0, 0.0]")
        if rot != [0.0, 0.0, 0.0]:
            self.ui.label("Mesh rotation " + str(rot) + " is not [0.0, 0.0, 0.0]")
        if scale != [0.0, 0.0, 0.0]:
            self.ui.label("Mesh scaling " + str(scale) + " is not [1.0, 1.0, 1.0]")
        self.ui.ops_button(rel_data_path="default.locrotscale", button_label="Fix (Move) Mesh Position")

    def draw_main_protein_mesh_panel(self):
        """
        Draw the main protein mesh panel.
        """

        # The location, rotation, and scaling are ok, so show normal UI
        self.ui.use_box_row("Load a Protein Trajectory")
        self.ui.object_property(property_name="pdb_filename")
        self.ui.new_row()

        self.ui.use_box_row("Simplify Trajectory")
        self.ui.object_property(property_name="frame_stride")
        self.ui.object_property(property_name="overall_pruning_stride")
        self.ui.new_row()

        self.ui.use_box_row("High-Detail Regions")

        # Go through and find the high-detail regions, list them.
        spheres = [obj for obj in bpy.data.objects if obj.name.startswith("Mineral_highres_sphere__")]
        for i, obj in enumerate(spheres[:10]):  # At most 10 displayed
            self.ui.ops_button(rel_data_path="select.sphere" + str(i), button_label="Sphere #" + str(i + 1) + " (Keep Every " + str(obj.sphere_pruning_stride) + " Atoms)")

        Messages.display_message("SELECT_SPHERE", self)
        self.ui.ops_button(rel_data_path="add.sphere", button_label="Create Region")

        # Create the button to start importing the MD trajectory.
        self.ui.use_box_row("Finalize")
        self.ui.label("WARNING: Loading simulation may take a bit.")
        Messages.display_message("TRAJ_FILENAME_DOESNT_EXIST", self)
        self.ui.ops_button(rel_data_path="load.traj", button_label="Load Trajectory")

        self.ui.new_row()

def menu_func(self, context):
    """
    Adds Mineral to Blender's menu system.

    :param bpy_types.Context context: The context.
    """

    self.layout.operator(Mineral.bl_idname)


class OBJECT_OT_LoadTrajButton(ButtonParentClass):
    """
    Button for loading in a trajectory.
    """

    bl_idname = "load.traj"
    bl_label = "Load Trajectory"

    def __init__(self):
        """
        Initialize the button object.
        """

        self.trajectory = None
        self.kdtree = None
        self.overall_pruning_stride = 1
        self.frame_stride = None

    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        global currently_loading_traj

        obj = context.object

        # Check to make sure filename exists
        if not os.path.exists(obj.pdb_filename):
            Messages.send_message("TRAJ_FILENAME_DOESNT_EXIST", "ERROR: Trajectory filename doesn't exist!")
        else:
            # Trajectory filename does exist, so load it...
            self.frame_stride = obj.frame_stride
            self.overall_pruning_stride = obj.overall_pruning_stride
            currently_loading_traj = True
            bpy.ops.process.trajectory('INVOKE_DEFAULT')

        return {'FINISHED'}

def geometric_center(obj):
    """
    Calculate the geometric center of an object.

    :param obj bpy_types.Object: The object to consider.

    :returns: a vector, the geometric center
    :rtype: :class:`Vector`
    """
    
    # See https://blender.stackexchange.com/questions/62040/get-center-of-geometry-of-an-object
    local_bbox_center = 0.125 * sum((Vector(b) for b in obj.bound_box), Vector())
    global_bbox_center = obj.matrix_world * local_bbox_center
    return global_bbox_center


class OBJECT_OT_AddSphereButton(ButtonParentClass):
    """
    Button for adding a high-detail sphere.
    """

    bl_idname = "add.sphere"
    bl_label = "Add Selection Sphere"

    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        obj = context.object

        # So dumb that blender throws an error if it's already in object
        # mode...
        try: bpy.ops.object.mode_set(mode='OBJECT')
        except: pass

        # Is the cursor near the mesh?
        cursor_loc = bpy.context.scene.cursor_location
        margin = 5.0
        global_bbox_center = geometric_center(obj)
        a_min = global_bbox_center - 0.5 * obj.dimensions
        a_max = global_bbox_center + 0.5 * obj.dimensions

        if (cursor_loc.x > a_min.x - margin and 
            cursor_loc.y > a_min.y - margin and 
            cursor_loc.z > a_min.z - margin and 
            cursor_loc.x < a_max.x + margin and 
            cursor_loc.y < a_max.y + margin and 
            cursor_loc.z < a_max.z + margin):

            # The 3D cursor is near the protein mesh. Add sphere at cursor
            # location.
            bpy.ops.mesh.primitive_uv_sphere_add(
                segments=16, 
                ring_count=16, 
                size=5.0, # Radius
                view_align=False, 
                enter_editmode=False, 
                location=cursor_loc, 
                rotation=(0.0, 0.0, 0.0)
            )
            
            # Store the new sphere in a variable.
            sphere = bpy.context.scene.objects.active

            # Pick sphere name, making sure not already used.
            sphere_name = "Mineral_highres_sphere__" + obj.name + "__" + str(0)
            i = 0
            while sphere_name in bpy.data.objects.keys():
                i = i + 1
                sphere_name = "Mineral_highres_sphere__" + obj.name + "__" + str(i)

            sphere.name = sphere_name

            # The sphere should be wireframe (to see into).
            bpy.ops.object.modifier_add(type='WIREFRAME')
            sphere.modifiers["Wireframe"].thickness = 0.2
        else:
            # The 3D cursor is not near the selected mesh, so throw an error...
            Messages.send_message("SELECT_SPHERE", "ERROR: Click on protein mesh to position 3D cursor!")
        return{'FINISHED'}


class OBJECT_OT_SphereDoneButton(ButtonParentClass):
    """
    Button to return to Protein Mesh panel once done with high-detail sphere.
    """

    bl_idname = "backto.protein"
    bl_label = "Back to Protein Mesh"

    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        # Figure out what the protein mesh is.
        obj = context.object
        protein_mesh_name = obj.name.split("__")[1]
        protein_mesh = bpy.data.objects[protein_mesh_name]

        # Make sure sphere pruning factor is less than whole protein pruning
        if protein_mesh.overall_pruning_stride < obj.sphere_pruning_stride:
            # It isn't, so throw an error.
            Messages.send_message("SPHERE_STRIDE_TOO_HIGH", "ERROR: Value too big. Set <= " + str(protein_mesh.overall_pruning_stride) + " (general value).")
        else:
            # It is, so switch back to the protein mesh.
            for obj in bpy.data.objects: obj.select = False

            bpy.context.scene.objects.active = protein_mesh
            bpy.context.scene.objects.active.select = True
        
        return{'FINISHED'}


class OBJECT_OT_DeleteSphereButton(ButtonParentClass):
    """
    Button to delete a sphere.
    """
    
    bl_idname = "delete.region"
    bl_label = "Delete Region"

    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        # Get the associated protein mesh.
        obj = context.object
        protein_mesh_name = obj.name.split("__")[1]
        protein_mesh = bpy.data.objects[protein_mesh_name]

        # Delete the sphere.
        for o in bpy.data.objects: o.select = False
        obj.select = True
        bpy.ops.object.delete() 

        # Switch back to the protein mesh.
        bpy.context.scene.objects.active = protein_mesh
        bpy.context.scene.objects.active.select = True

        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButtonParent(ButtonParentClass):
    """
    Button to select and edit an existing sphere. This is a parent all other
    select-sphere buttons inherit. 
    """

    bl_idname = ""
    bl_label = ""

    def switch_to_obj(self, index):
        """
        Switch to a given high-detail sphere.

        :param int index: The index of the sphere.
        """
        
        # Get the sphere
        spheres = [obj for obj in bpy.data.objects if obj.name.startswith("Mineral_highres_sphere__")]
        sphere = spheres[index]

        # Make that sphere selected and active.
        for obj in bpy.data.objects: obj.select = False
        bpy.context.scene.objects.active = bpy.data.objects[sphere.name]
        bpy.context.scene.objects.active.select = True


class OBJECT_OT_SelectExistingSphereButton0(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere0"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(0)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton1(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere1"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(1)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton2(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere2"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(2)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton3(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere3"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(3)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton4(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere4"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(4)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton5(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere5"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(5)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton6(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere6"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(6)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton7(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere7"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(7)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton8(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere8"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(8)
        return{'FINISHED'}


class OBJECT_OT_SelectExistingSphereButton9(OBJECT_OT_SelectExistingSphereButtonParent):
    """
    Button to select and edit an existing sphere.
    """

    bl_idname = "select.sphere9"
    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        self.switch_to_obj(9)
        return{'FINISHED'}


class OBJECT_OT_StartOver(ButtonParentClass):
    """
    Button to start over entirely. Deletes all animations and high-detail
    spheres.
    """

    bl_idname = "start.over"
    bl_label = "Start Over"

    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        # Delete anything object that starts with the plugin's name.
        bpy.ops.object.select_all(action='DESELECT')
        for obj in bpy.data.objects:
            if obj.name.startswith("Mineral_"):
                obj.select = True
                bpy.ops.object.delete()

        return{'FINISHED'}


class OBJECT_OT_RemoveAnimations(ButtonParentClass):
    """
    Button to start over in part. Removes animations, but leaves high-detail
    spheres intact.
    """

    bl_idname = "remove.animations"
    bl_label = "Remove Animations"

    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        # Delete anything that starts with the plugin's name and is not a
        # high-detail sphere.
        bpy.ops.object.select_all(action='DESELECT')
        for obj in bpy.data.objects:
            if (
                obj.name.startswith("Mineral_") and 
                not obj.name.startswith("Mineral_highres_sphere__")
            ):
                obj.select = True
                bpy.ops.object.delete()

        return{'FINISHED'}


class OBJECT_OT_DefaultLocRotScaleButton(ButtonParentClass):
    """
    Button for setting a protein mesh's location, rotation, and scaling
    vectors.
    """

    bl_idname = "default.locrotscale"
    bl_label = "Fix (Move) Mesh Position"

    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        obj = context.object

        # So dumb that blender throws an error if it's already in object mode...
        try: bpy.ops.object.mode_set(mode='OBJECT')
        except: pass

        # Set the object's location.
        bpy.context.scene.objects.active.location.x = 0.0
        bpy.context.scene.objects.active.location.y = 0.0
        bpy.context.scene.objects.active.location.z = 0.0

        # Set the object's rotation.
        bpy.context.scene.objects.active.rotation_euler.x = 0.0
        bpy.context.scene.objects.active.rotation_euler.y = 0.0
        bpy.context.scene.objects.active.rotation_euler.z = 0.0

        # Set the object's scale.
        bpy.context.scene.objects.active.scale.x = 1.0
        bpy.context.scene.objects.active.scale.y = 1.0
        bpy.context.scene.objects.active.scale.z = 1.0
        
        # Zoom in on the object.
        for obj in bpy.data.objects: obj.select = False
        bpy.context.scene.objects.active.select = True
        bpy.ops.view3d.view_selected(use_all_regions=False)

        return{'FINISHED'}


class OBJECT_OT_MainMenuButton(ButtonParentClass):
    """
    Button to return to the main menu.
    """

    bl_idname = "main.menu"
    bl_label = "Return to Main Menu"

    def execute(self, context):
        """
        Runs when button pressed.

        :param bpy_types.Context context: The context.
        """

        obj = context.object

        # So dumb that blender throws an error if it's already in object
        # mode...
        try: bpy.ops.object.mode_set(mode='OBJECT')
        except: pass

        for obj in bpy.data.objects:
            obj.select = False
        
        return {'FINISHED'}

# Store keymaps here to access after registration.
addon_keymaps = []
classes_used = [
    OBJECT_OT_LoadTrajButton,
    OBJECT_OT_AddSphereButton,
    OBJECT_OT_DefaultLocRotScaleButton,
    OBJECT_OT_SphereDoneButton,
    OBJECT_OT_SelectExistingSphereButton0,
    OBJECT_OT_SelectExistingSphereButton1,
    OBJECT_OT_SelectExistingSphereButton2,
    OBJECT_OT_SelectExistingSphereButton3,
    OBJECT_OT_SelectExistingSphereButton4,
    OBJECT_OT_SelectExistingSphereButton5,
    OBJECT_OT_SelectExistingSphereButton6,
    OBJECT_OT_SelectExistingSphereButton7,
    OBJECT_OT_SelectExistingSphereButton8,
    OBJECT_OT_SelectExistingSphereButton9,
    OBJECT_OT_DeleteSphereButton,
    OBJECT_OT_StartOver,
    OBJECT_OT_RemoveAnimations,
    ProcessTrajectory,
    OBJECT_OT_MainMenuButton
]

##### Registration functions #####
def register():
    """
    Registers this addon.
    """
    Mineral.start()
    bpy.utils.register_class(Mineral)
    bpy.types.VIEW3D_MT_object.append(menu_func)

    global classes_used
    for c in classes_used:
        bpy.utils.register_class(c)

def unregister():
    """
    Good practice to make it possible to unregister addons.
    """

    bpy.utils.unregister_class(Mineral)
    bpy.types.VIEW3D_MT_object.remove(menu_func)

    global classes_used
    for c in classes_used:
        bpy.utils.unregister_class(c)

if __name__ == "__main__":
    """
    Start the addon!
    """

    register()
