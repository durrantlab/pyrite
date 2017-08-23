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

import bpy
import time

class BackgroundJobParentClass(bpy.types.Operator):
    """
    A class for running jobs in the background.
    """

    bl_idname = "object.modal_operator"
    bl_label = ""

    first_run = True

    def modal(self, context, event):  # The main loop
        """
        The main loop. This runs periodically to drive the job. So the job
        needs to be divided into steps.

        :param bpy_types.Context context: The context.

        :param bpy.types.Event event: The event.
        """

        # In case the user cancels the job.
        if event.type in {'ESC'}:
            self.job_cancelled()
            return {'CANCELLED'}
        
        # First time, setup stuff.
        if self.first_run == True:
            self.first_run = False
            self.setup(context, event)

        # Run a single test.
        response = self.run_step(context, event)
        if response is not None:
            return response

        # return {'RUNNING_MODAL'}  # If you only want it to run when event fires.
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        """
        I believe this runs when the class is loaded, not when first called.

        :param bpy_types.Context context: The context.

        :param bpy.types.Event event: The event.
        """

        if context.object:
            self.timer = context.window_manager.event_timer_add(0.01, context.window)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "No active object, could not finish")
            return {'CANCELLED'}

    def job_cancelled(self):
        """
        The user cancels the job. Overwritten by children.
        """

        pass

    def setup(self, context, event):
        """
        Set up variables needed to process the trajectory. Overwritten by children.

        :param bpy_types.Context context: The context.

        :param bpy.types.Event event: The event.
        """
        
        pass

    def run_step(self, context, event):
        """
        Run a single step of the background job. You need to do it in steps so
        you can periodically return control to the UI, making it seem like a
        background job. Overwritten by children.

        :param bpy_types.Context context: The context. 

        :param bpy.types.Event event: The event.
        """

        pass

