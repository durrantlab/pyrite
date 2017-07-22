import bpy
 
class OBJECT_PT_pingpong(bpy.types.Panel):
    bl_label = "Ping Pong"
    bl_space_type = "PROPERTIES"
    bl_region_type = "TOOLS"
    bl_context = "object"
 
    is_left = True
 
    def draw_header(self, context):
        layout = self.layout
        layout.label(text="", icon="PHYSICS")
 
    def draw(self, context):
        layout = self.layout
 
        row = layout.row()
        split = row.split(percentage=0.5)
        col_left = split.column()
        col_right = split.column()
 
        if self.is_left:
            col_left.operator("object.pingpong", text="Ping")
        else:
            col_right.operator("object.pingpong", text="Pong")
 
 
class OBJECT_OT_pingpong(bpy.types.Operator):
    bl_label = "Ping Pong Operator"
    bl_idname = "object.pingpong"
    bl_description = "Move the ball"
 
    def execute(self, context):
        OBJECT_PT_pingpong.is_left = not OBJECT_PT_pingpong.is_left
        self.report({'INFO'}, "Moving the ball")
        return {'FINISHED'}
 
def register():
    bpy.utils.register_module(__name__)
 
def unregister():
    bpy.utils.unregister_module(__name__)
 
if __name__ == "__main__":
    register()