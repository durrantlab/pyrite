# DurBlend is a library for easily creating blender plugins. 
# Copyright (C) 2017  Jacob D. Durrant
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

import textwrap
import bpy

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
        if label_txt != "":
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

    def multiline_label(self, txt="moose dog face", col_width=15):
        """Never got this to work..."""
        """
        # Get width of current area.
        areas = bpy.context.screen.areas
        for i, a in enumerate(areas):
            if (a.type == "VIEW_3D"):
                aID = i
            reg = a.regions
            for ir, r in enumerate(reg):
                if (r.type == "TOOL_PROPS"):
                    rID = ir
                    print(aID, rID, "MOO")
        if ((aID >= 0) and (rID >= 0)):
            width = areas[aID].regions[rID].width

            # Determine the number of letters per line. For some reason, this
            # doesn't seem to scale linearly. Also, I suspect this will depend
            # on things like font size... oh well.
            if width < 425:
                pixels_per_letter = 20
            elif width < 525:
                pixels_per_letter = 17
            elif width < 925:
                pixels_per_letter = 15
            else:
                pixels_per_letter = 13
            col_width = int(width / pixels_per_letter)
        else:
            # Couldn't get panel width... not sure what's going on... so just
            # show label.
            self.label(txt)
            return

        for line_index, line in enumerate(textwrap.wrap(txt, col_width)):
            # num_spaces_to_add = col_width - len(line)
            # half_num_spaces_to_add = 0.5 * num_spaces_to_add
            # add_to_left = " " * int(half_num_spaces_to_add)
            # if line_index % 2 == 1:
            #     add_to_left = [" ", ""][int(0.5 * (line_index - 1)) % 2] + add_to_left

            # line = (add_to_left + line).ljust(col_width)

            row = self.new_row()
            # print(dir(row))
            # areas=bpy.context.screen.areas
            # print(areas)
            # p1 = row.split(percentage=10)
            row.label(text=line)
            row.scale_y = 0.55
            # row.alignment='CENTER'
        """

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
        row.operator(rel_data_path, text=button_label)

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
