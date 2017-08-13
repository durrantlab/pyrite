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

try: import cStringIO as StringIO
except: from io import StringIO

from . import scoria
import numpy
import mathutils
import bpy
from bpy import context
from mathutils import Vector
import os
import tempfile
import glob
import shutil

from .DurBlend import Properties
from .DurBlend import UI
from .DurBlend import PanelParentClass
from .DurBlend import ButtonParentClass
from .DurBlend import Messages
from .DurBlend import BackgroundJobParentClass

plugin_name = "Mineral"

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

        global plugin_name
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
            if not obj_to_use.name.startswith(plugin_name + "_highres_sphere__")
            else obj_to_use.name.split("__")[1]
        )

        if obj_to_use.name.startswith(plugin_name + "_highres_sphere__"):
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
            if obj.name.startswith(plugin_name + "_"):
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
        self.ui.label("If you use " + plugin_name + ", please cite:")
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
        spheres = [obj for obj in bpy.data.objects if obj.name.startswith(plugin_name + "_highres_sphere__")]
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

        global plugin_name

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
        global plugin_name
        
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
            sphere_name = plugin_name + "_highres_sphere__" + obj.name + "__" + str(0)
            i = 0
            while sphere_name in bpy.data.objects.keys():
                i = i + 1
                sphere_name = plugin_name + "_highres_sphere__" + obj.name + "__" + str(i)

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
        
        global plugin_name

        # Get the sphere
        spheres = [obj for obj in bpy.data.objects if obj.name.startswith(plugin_name + "_highres_sphere__")]
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

        global plugin_name

        # Delete anything object that starts with the plugin's name.
        bpy.ops.object.select_all(action='DESELECT')
        for obj in bpy.data.objects:
            if obj.name.startswith(plugin_name + "_"):
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

        global plugin_name

        # Delete anything that starts with the plugin's name and is not a
        # high-detail sphere.
        bpy.ops.object.select_all(action='DESELECT')
        for obj in bpy.data.objects:
            if (
                obj.name.startswith(plugin_name + "_") and 
                not obj.name.startswith(plugin_name + "_highres_sphere__")
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


class ProcessTrajectory(BackgroundJobParentClass):
    """
    A class to load and process a molecular dynamics trajectory.
    """
    
    bl_idname = "process.trajectory"

    def setup(self, context, event):
        """
        Set up variables needed to process the trajectory.

        :param bpy_types.Context context: The context.

        :param bpy.types.Event event: The event.
        """
        
        self.overall_pruning_stride = context.object.overall_pruning_stride
        self.current_step = "START"
        self.current_frame = None
        self.pruning_spheres = []
        self.selection_atoms_to_keep_with_offset = {}
        self.selection_atoms_to_keep = []
        self.frame_index = 0
        self.frame_stride = context.object.frame_stride
        self.pdb_filename = context.object.pdb_filename
        self.frames = self.get_frames()
        self.guide_empties_location = None
        self.protein_obj = bpy.context.scene.objects.active

    def run_step(self, context, event):
        """
        Run a single step of the background job. You need to do it in steps so
        you can periodically return control to the UI, making it seem like a
        background job.

        :param bpy_types.Context context: The context. 

        :param bpy.types.Event event: The event.
        """
        
        if self.current_step == "START":
            # On the first step, set the current frame and post a message.
            self.current_frame = next(self.frames)
            Messages.send_message(
                "LOAD_TRAJ_PROGRESS", 
                "Identifying which atoms to keep..."
            )
            self.current_step = "ID_ATOMS"
        elif self.current_step == "ID_ATOMS":
            # Second step. Loading the first frame.
            self.get_pruned_indecies()
            self.add_empties_at_atom_points()
            Messages.send_message("LOAD_TRAJ_PROGRESS", "Loading frame 1...")
            self.current_step = "POSITION_EMPTIES"
        if self.current_step == "LOAD_FRAME":
            # First half if third step. Loading additional frames...
            try:
                self.current_frame = next(self.frames)
                Messages.send_message(
                    "LOAD_TRAJ_PROGRESS", 
                    "Loading frame " + str(self.frame_index + 1) + "..."
                )
                self.current_step = "POSITION_EMPTIES"
                return None
            except StopIteration:
                # Here you've reached the last frame, so proceed to the next
                # step.
                self.position_empties_at_atom_locs(position_all=True)
                Messages.send_message(
                    "LOAD_TRAJ_PROGRESS", 
                    "Setting up armature and bones..."
                )
                self.current_step = "ARMATURE"
        elif self.current_step == "POSITION_EMPTIES":
            # Second half of third step. Position empties at the location of
            # the atoms.
            self.position_empties_at_atom_locs()
            self.frame_index = self.frame_index + 1
            self.current_step = "LOAD_FRAME"
        elif self.current_step == "ARMATURE":
            # Everything loaded. Now set up armature and bones.
            self.armature_and_bones()
            return {'CANCELLED'}
    
    def job_cancelled(self):
        """
        The user cancels the job.
        """

        global currently_loading_traj
        bpy.ops.remove.animations('INVOKE_DEFAULT')
        currently_loading_traj = False

    def get_frames(self):
        """
        Get the next frame from the trajectory. Because this is a generator,
        we only load in one frame at a time. That's good for memory.
        """

        current_frame = 0
        pdb_file = open(self.pdb_filename, 'r')
        current_pdb_frame_lines = []

        while True:
            line = pdb_file.readline()
            current_pdb_frame_lines.append(line)
            if line.strip() in ["END", "ENDMDL"]:
                print("Frame", current_frame)

                mol = scoria.Molecule()
                mol.load_pdb_into_using_file_object(
                    StringIO("".join(current_pdb_frame_lines)),
                    False, False, False, False
                )
                yield mol   

                # Start over with next frame
                current_pdb_frame_lines = []
                current_frame += 1

            if line == "":  # EOF
                break

    def add_empties_at_atom_points(self):
        """
        Add enough empty objects for the number of atoms that will be loaded.
        """

        # Make an empty that all the other empties will be parented to. This
        # is just to keep things organized.
        guide_empties = bpy.ops.object.empty_add(type='PLAIN_AXES', radius=1)
        guide_empties = bpy.context.object
        guide_empties.name = plugin_name + "_GuideEmpties"
        self.guide_empties_location = numpy.mean(
            self.current_frame.get_coordinates(), 
            axis=0
        )
        guide_empties.location = self.guide_empties_location
        
        # Add enough empties to match the number of bones. 
        for index in self.selection_atoms_to_keep:
            empty = bpy.ops.object.empty_add(type='PLAIN_AXES', radius=1)
            empty = bpy.context.object
            empty.name = plugin_name + "_empty" + str(index)
            empty.parent = guide_empties

    def get_pruned_indecies(self):
        """
        Applies pruning spheres to the existing protein. Determines the
        indecies of the atoms to keep.
        """

        # The key is to use the smallest pruning stride possible for a given
        # point.

        global plugin_name

        # Add the pruning spheres
        # Add oveall pruning sphere
        self.pruning_spheres.append((self.overall_pruning_stride, 0.0, 0.0, 0.0, 1e50))

        # Add individual spheres
        for sphere in [obj for obj in bpy.data.objects if obj.name.startswith(plugin_name + "_highres_sphere__")]:
            x, y, z = list(sphere.location)
            r = 5.0 * sphere.scale.x  # Because radius is 5 at creation
            stride = sphere.sphere_pruning_stride
            self.pruning_spheres.append((stride, x, y, z, r))

        # Make a kd tree if needed
        coors = self.current_frame.get_coordinates()  # So kdtree calculated on
                                                      # coordinates of first frame only.
        # Create kd-tree containing atom/bone coordinates
        kdtree = mathutils.kdtree.KDTree(len(coors))
        bone_list = []

        # Add coordinates to tree
        for i, c in enumerate(coors):
            kdtree.insert(c, i)
        kdtree.balance()

        # Make sure the pruning spheres are ordered by the stride, from
        # smallest stride to greatest.
        self.pruning_spheres.sort()

        # Go through each sphere and apply a mask, where true means the
        # coordinate is in the given sphere, and false means it isn't.
        indices_in_previous_spheres = set([])
        total_indices_to_keep = set([])
        for sphere in self.pruning_spheres:
            # Get the coordinate indices that are in the sphere
            atom_stride, center_x, center_y, center_z, radius = sphere
            co_find = (center_x, center_y, center_z)
            coors_in_sphere = kdtree.find_range(co_find, radius)
            # coors_in_sphere looks like this:
            # (Vector((-26.695999145507812, -92.79199981689453, 389.52801513671875)), 777, 10.07503890991211)

            # Just keep the indices in this sphere
            indices_in_this_sphere = set([])
            for coor in coors_in_sphere:
                indices_in_this_sphere.add(int(coor[1]))
            
            # Remove from the ones in this sphere ones that were in previous spheres
            indices_to_keep = indices_in_this_sphere - indices_in_previous_spheres

            # Update indices_in_previous_spheres
            indices_in_previous_spheres = indices_in_previous_spheres.union(indices_to_keep)

            # Now only keep every few of the indices_to_keep
            indices_to_keep_sparce = numpy.array(list(indices_to_keep))[::atom_stride]

            # Update all the totals
            total_indices_to_keep = total_indices_to_keep.union(set(indices_to_keep_sparce))

        # Keep track of the atomic indecies to keep.
        total_indices_to_keep = list(total_indices_to_keep)
        total_indices_to_keep.sort()
        total_indices_to_keep = numpy.array(total_indices_to_keep)
        self.selection_atoms_to_keep = total_indices_to_keep
        for i in range(self.frame_stride):
            self.selection_atoms_to_keep_with_offset[i] = total_indices_to_keep[i::self.frame_stride]

    def position_empties_at_atom_locs(self, position_all=False):
        """
        Positions the empties at the locations of the retained atoms.

        :param boolean position_all: Whether or not to position all empties.
                       Use True at first and last frame. Defaults to False.
        """

        try:  
            # So dumb that blender throws an error if it's already in object
            # mode...
            bpy.ops.object.mode_set(mode='OBJECT')
        except:
            pass

        global plugin_name

        # Sets next frame to add
        bpy.context.scene.frame_set(self.frame_index)

        if self.frame_index == 0 or position_all == True:
            # Set all of them.
            sel = self.selection_atoms_to_keep
        else:
            # Set just the ones appropriate for this frame (using offset).
            sel = self.selection_atoms_to_keep_with_offset[
                self.frame_index % self.frame_stride
            ]
        
        # Get the coordinates of the appropriate atoms.
        coors = self.current_frame.get_coordinates()[sel]

        # Go through and set the appropriate empty object to the associated
        # coordinate.
        for coor_index, coor in enumerate(coors):
            si = sel[coor_index]
            empty = bpy.data.objects[plugin_name + "_empty" + str(si)]
            empty.location = coor - self.guide_empties_location
            empty.keyframe_insert(data_path='location')

    def armature_and_bones(self):
        """
        Setup the armature and bones.
        """

        global plugin_name
        global currently_loading_traj

        # Creating armature
        bpy.ops.object.add(type='ARMATURE', enter_editmode=True)
        object = bpy.context.object
        object.name = plugin_name + '_Armature'
        armature = object.data
        armature.name = 'Frame'

        # Add bones
        for index in self.selection_atoms_to_keep:
            bone_name = plugin_name + '_bone' + str(index)
            bone = armature.edit_bones.new(bone_name)
            bone.head = (0, 0, 0)
            bone.tail = (0, 0, 2)
            bone.envelope_weight = 1.0  # Needed for envelope-based mesh
                                        # vertex weighting.
            bone.envelope_distance = 5.0

        # Now constrain them to the locations of the empty objects (which are
        # animated).
        bpy.ops.object.mode_set(mode='POSE')
        armature = bpy.data.objects[plugin_name + "_Armature"]

        for index in self.selection_atoms_to_keep:
            bone = armature.pose.bones[plugin_name + '_bone' + str(index)]
            constraint = bone.constraints.new(type="COPY_LOCATION")
            constraint.target = bpy.data.objects[plugin_name + "_empty" + str(index)]

        # Now make sure the pose at frame 0 is set as the rest pose (so you
        # can do automatic weights later...)
        bpy.context.scene.frame_set(0)
        bpy.ops.pose.armature_apply()

        # Parent the protein mesh to that armature
        for obj in bpy.data.objects: obj.select = False
        self.protein_obj.select = True
        armature.select = True
        bpy.context.scene.objects.active = armature
        bpy.ops.object.parent_set(type="ARMATURE_AUTO")

        # Zoom in on armature
        try: bpy.ops.object.mode_set(mode='OBJECT')
        except: pass
        
        bpy.context.scene.objects.active = bpy.data.objects[plugin_name + '_Armature']
        bpy.ops.view3d.view_selected(use_all_regions=False)

        # No longer running
        currently_loading_traj = False


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
