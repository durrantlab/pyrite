Pyrite 
=======

"Pyrite" is a program for importing a molecular dynamics trajectory into
Blender, to take advantage of Blender's advanced rendering features. A copy of
the plugin is available at
[http://durrantlab.com/pyrite/](http://durrantlab.com/pyrite/), released under
the terms of the GNU General Public License Version 3.

If you use Pyrite, please cite our paper: [CITATION HERE]()


Tutorial
========

1. Create a Protein Mesh
------------------------

Free programs such as VMD, PyMOL, and Chimera export molecular representations
in Blender-compatible formats (e.g., OBJ, WRL, X3D, etc). Exported mesh files
use camera coordinates rather than the world coordinates of the original
model. Before exporting, all transformation matrices must be set to identity.
In VMD, this simple TCL script, adapted from code provided by John Stone
(VMD’s primary developer), sets the coordinate system: 

```
set m {
  {{1.0 0.0 0.0 0.0}
   {0.0 1.0 0.0 0.0}
   {0.0 0.0 1.0 0.0}
   {0.0 0.0 0.0 1.0}}
  {{1.0 0.0 0.0 0.0}
   {0.0 1.0 0.0 0.0}
   {0.0 0.0 1.0 0.0}
   {0.0 0.0 0.0 1.0}}
  {{1.0 0.0 0.0 0.0}
   {0.0 1.0 0.0 0.0}
   {0.0 0.0 1.0 0.0}
   {0.0 0.0 0.0 1.0}}
  {{1.0 0.0 0.0 0.0}
   {0.0 1.0 0.0 0.0}
   {0.0 0.0 1.0 0.0}
   {0.0 0.0 0.0 1.0}}
}

for {set i 0} {$i < [molinfo num]} {incr i} {
  molinfo ${i} set {center_matrix rotate_matrix scale_matrix global_matrix} $m
}
```

In PyMOL, the same can be accomplished with this Python script, which moves
and rotates the camera to the appropriate location:

```
cmd.set_view([ 1.0,   0.0,   0.0,
               0.0,   1.0,   0.0,
               0.0,   0.0,   1.0,
               0.0,   0.0,   0.0,
               0.0,   0.0,   0.0,
              40.0, 100.0, -20.0])
```

After setting up the scene in your program of choice (e.g., using ribbon or
surface representation for the proteins), export the first frame as an OBJ
file.

2. Import Protein Mesh into Blender
-----------------------------------

Import your Wavefront OBJ file using Blender's menu: ```File > Import >
Wavefront (.obj)```.  After importing the mesh, prepare it for animation. For
example, using Blender's “Remove Doubles” command is often critical. Available
in Edit Mode, this command merges duplicate vertices into one to ensure that
all mesh faces are connected.

Blender itself sometimes imports meshes using the incorrect coordinate system.
Setting the object position and rotation vectors to ```(0.0, 0.0, 0.0)```, and
the scaling vector to ```(1.0, 1.0, 1.0)```, ensures that the imported object
and trajectory use the same system. Pressing the "Auto-Fix" button in Pyrite's
Tool Shelf Panel will automatically set these vectors. You can also set the
values in Blender’s Object Properties Panel.

3. Specify the Trajectory File
------------------------------

If an object has the appropriate transformation vectors, the Pyrite tab
presents options for simulation import. Use the first option to specify the
location of the MD trajectory file, in multi-frame PDB format. Trajectories
saved in other formats (e.g., the binary DCD format) can be converted to PDB
using [VMD](http://www.ks.uiuc.edu/Research/vmd/) or
[catdcd](http://www.ks.uiuc.edu/Development/MDTools/catdcd/). 

4. Set the Trajectory-Simplification Options
--------------------------------------------

Loading the position of every atom across the entire MD trajectory is often
too memory and CPU intensive. Instead, Pyrite coarse grains the simulation
across time and space by discarding some frames and atoms. Select the temporal
and spatial resolution by indicating how often to keep a frame or atom (e.g.,
only every nth frame and every mth atom). 

5. Create High-Detail Regions
-----------------------------

You may wish to more accurately represent the MD-captured motions of some
regions (e.g., active sites):

1. Add a sphere to Blender's viewport by positioning the 3D cursor and
   clicking the “Create Region” button in the “High-Detail Regions” subpanel. 
2. Adjust the position and scaling until the sphere encompasses the region of
   interest. 
3. Use the Pyrite panel to specify how often to keep atoms within the sphere.
   By selecting a lower skip/stride value, fewer atoms within the
   sphere-defined region are discarded. The captured MD motions within this
   region will have higher spatial resolution.
4. After positioning and scaling the sphere, press the “Back to Protein Mesh”
   button to save. Alternatively, the “Delete Region” button will remove the
   current sphere. 

You can add several spheres to the scene if multiple distinct regions require
higher spatial resolution.

6. Load the Trajectory
----------------------

After saving the high-detail regions, press the “Load Trajectory” button to
import trajectory data. Empty objects populate the viewport at the locations
of the retained atoms. 

