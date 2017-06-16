import bpy
from bpy.props import *

##### To test this plugin, paste it into Blender's text editor and click "Run
#Script." #####

##### Plugin information #####
bl_info = {
    "name" : "Tetrahedron Maker",
    "author" : "J Random Hacker <jrhacker@example.com>",
    "version" : (1, 0, 0),
    "blender" : (2, 5, 7),
    "location" : "View 3D > Edit Mode > Tool Shelf",
    "description" : "Generate a tetrahedron mesh",
    "warning" : "",
    "wiki_url" : "",
    "tracker_url" : "",
    "category" : "Add Mesh",
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
            name = txt,
            min = min, max = max,
            default = default,
            description = "An integer between " + str(min) + " and " + str(max),
            update = update
        )
    
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
            name = txt,
            min = min, max = max,
            default = default,
            description = "A float between " + str(min) + " and " + str(max),
            update = update
        )
    
    def boolProp(self, txt, default = True, update = nothing):
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
            name = txt,
            default = default,
            description = "True or false",
            update = update
        )
    
    def strProp(self, txt, default="", update=nothing):
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
            default = default,
            description = "Text",
            update=update
        )
    
    def enumProp(self, txt, items=[("moose", "Moose", ""), ("dog", "Dog", "")], update = nothing):
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
            name = txt,
            #default = items[0],
            description = "Select Option",
            update = update,
            items = items
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
        your scene, you could use this button with the rel_daata_path set to
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
        the rel_daata_path set to "object.select_all" won't do. You
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

        self.setup_properties()  # Can you call a class def like this? If you get an error, you can't.

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


##### An Example Panel #####
class GoatsAreGreat(PanelParentClass):
    bl_label = "Goats Are Great"
    bl_idname = "MATERIALS_PT_greatgoats"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    
    @classmethod
    def setup_properties(self):
        """
        Define all the scene and object properties for your panel. Every Panel
        class must have this function!
        """

        # Set up scene and object properties.
        bpy.types.Object.goat_name = self.prop_funcs.strProp("Name", "Goaty McGoatFace")
        bpy.types.Object.is_goat_great = self.prop_funcs.boolProp("Goat is great", False, update = self.check_goat_opinion)
        bpy.types.Object.rating = self.prop_funcs.floatProp("Probability goat is pink", min = 0.0, max = 1.0, default = 0.5)

        bpy.types.Object.num_legs = self.prop_funcs.intProp("# legs", min = 0, max = 4, default = 4)
        bpy.types.Object.goat_type = self.prop_funcs.enumProp("Goat Type", items=[("fainting", "Fainting", "F"), ("billy", "Billy", "B"), ("mountain", "Mountain", "M")])
        
        # You can also store variables internally.
        bpy.types.Object.user_opinion = self.prop_funcs.strProp("", "I am eager to hear your opinion!")

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
        self.ui.label("Goat (Object Name: " + self.obj.name + ")")
        
        self.ui.use_box_row("Goat Identity")
        self.ui.object_property(property_name="goat_name")

        self.ui.use_box_row("Goat Greatness?")
        self.ui.object_property(property_name="is_goat_great")
        self.ui.object_property(property_name="rating")
        self.ui.label(context.object["user_opinion"])
        
        self.ui.use_box_row("About this Goat")
        self.ui.object_property(property_name="num_legs")
        self.ui.object_property(property_name="goat_type")

        self.ui.use_box_row("Goats Love Buttons")
        self.ui.ops_button(rel_data_path="object.modifier_add", button_label="Add a Modifier (ops_button)")
        self.ui.ops_action_button(rel_data_path="object.select_all", button_label="Invert Selection (ops_action_button)", action="INVERT")
        
    def check_goat_opinion(self, context):
        """
        This shows how you can do something extra when a value changes, using
        a callback function. It validates the users opinion when they click
        the Goats are Great checkbox.

        :param ??? context: The context of the currently selected object.
        """

        obj = context.object
        if obj['is_goat_great'] == False:
            msg = "False! All goats, great in their own way."
            print(msg)
            obj['user_opinion'] = msg
        else:
            msg = "You get the greatness of goats!"
            print(msg)
            obj['user_opinion'] = msg

##### Custom operators #####
# Can be called from the parameter
# Need to work on this more...
# class runScript(bpy.types.Operator):
#     bl_idname = "object.run_script"
#     bl_label = "Invokes a Script"
#     bl_options = {"UNDO"}

#     @classmethod
#     def poll(cls, context):
#         return context.active_object is not None

#     def execute(self, context):
#         print('Yay!')

#         return {'FINISHED'}


##### Registration functions #####
def register():
    """
    Registers this addon. MAKE SURE TO USE YOUR CLASS BELOW, NOT
    GoatsAreGreat!!!
    """

    GoatsAreGreat.start()
    bpy.utils.register_class(GoatsAreGreat)

def unregister():
    """
    Good practice to make it possible to unregister addons.
    """

    bpy.utils.unregister_class(__name__)

if __name__ == "__main__":
    """
    Start the addon!
    """

    register()