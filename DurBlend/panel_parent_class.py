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

import bpy
from .properties import Properties
from .ui import UI

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
        """
        Set certain variables.

        :param bpy_types.Context context: The context.
        """

        self.obj = context.object
        self.scene = bpy.context.scene
        self.ui.parent = self
