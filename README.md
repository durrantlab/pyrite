Nivedita's project 
==================

Note that the pymolecule version here is not necessarily the latest version.

Testing Plugin in Python
========================
To test the plugin in python, open up a Blender python console and add the
plugin directory to the python path:

`import sys
sys.path.append("/Users/jdurrant/Documents/Work/durrant_git/nivedita/blender_addon")`

Then import the module and register the plugin:

`import testing_addon
testing_addon.register()
`

If you make changes to the plugin, you need to unregister, reload the module, and reregister.:

`testing_addon.unregister()
import imp
imp.reload(testing_addon)
testing_addon.register()
`


