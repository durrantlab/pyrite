import bpy
import time

class BackgroundJobParentClass(bpy.types.Operator):
    bl_idname = "object.modal_operator"
    bl_label = ""

    first_run = True

    def modal(self, context, event):  # The main loop
        if event.type in {'ESC'}:
            print("Cancel")
            return {'CANCELLED'}
        
        # First time, setup stuff.
        if self.first_run == True:
            self.first_run = False
            self.setup(context, event)

        response = self.run_step(context, event)
        if response is not None:
            return response


        # return {'RUNNING_MODAL'}  # If you only want it to run when event fires
        return {'PASS_THROUGH'}

    def invoke(self, context, event):  # When loaded, I think. Not when first called.
        if context.object:
            # self.source_obj = context.object
            # self.setup()
            self.timer = context.window_manager.event_timer_add(0.01, context.window)
            # self.timer = context.window_manager.event_timer_add(1.0, context.window)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "No active object, could not finish")
            return {'CANCELLED'}

    # def start(self, **kwargs):  # My start function
    #     bpy.ops.object.modal_operator('INVOKE_DEFAULT')
    #     self.setup(kwargs)
    
    def setup(self, context, event):
        # Overwritten by children
        pass

    def run_step(self, context, event):
        # Overwritten by children
        pass

