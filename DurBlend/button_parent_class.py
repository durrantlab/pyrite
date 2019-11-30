# DurBlend is a library for easily creating blender plugins.
# Copyright (C) 2018  Jacob D. Durrant
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


class ButtonParentClass(bpy.types.Operator):
    """
    The parent class doesn't do anything for now. But in the future I might
    want to make certain functions common to all buttons. So all buttons
    inherit this one.
    """

    def _nothing(self): pass
