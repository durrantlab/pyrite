Pyrite 
=======

"Pyrite" is a program for importing a molecular dynamics trajectory into
Blender, to take advantage of Blender's advanced rendering features. The basic
steps are:

1. Import a protein mesh (created using VMD's Render feature). There are
   sample files in `./blender_addon/trajectory_samples/`
2. Specify the trajectory file in multi-frame PDB format (again, see
   `./blender_addon/trajectory_samples/`)
3. General coarse graining: Specify how often frames should be dropped, and
   how many atoms should be retained.
4. Use spheres to specify high-detail regions. The atom stride should be
   smaller in these regions (i.e., finer coarse grained).
5. Load in the specified atoms ("bones")/frames. Interpolate between frames to
   animate.
6. Skin these atoms/bones with the protein mesh.

What's Not Yet Implemented
==========================

I believe step 4 above isn't yet fully implemented. Also, the UI is
incomplete.

Testing the Python Plugin in Blender
====================================
To test the plugin in python, open up a Blender python console and add the
plugin directory to the python path:

```
import sys
sys.path.append("/Users/jdurrant/Documents/Work/durrant_git/nivedita/blender_addon")
```

Then import the module and register the plugin:

```
import pyrite
pyrite.register()
```

If you make changes to the plugin, you need to reload the module and
reregister.:

```
import imp
imp.reload(pyrite)
pyrite.register()
```

Misc
====

Note that the scoria version here is not necessarily the latest version.

