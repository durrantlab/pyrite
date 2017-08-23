# DurBlend is a library for easily creating blender plugins.
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

from bpy.props import *

def nothing(self, context):
    """
    This is a function that does nothing.
    """
    return

##### Setup scene and object variables #####
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

