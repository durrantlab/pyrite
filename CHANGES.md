Changes
=======

1.1.1
-----

* Fixed a bug that prevented Pyrite from loading some multi-frame PDB
  trajectories. Programs such as PyMOL separate frames using ENDMDL. Previous
  versions of Pyrite worked best with VMD-formatted files, which separate
  frames using END only.

1.1.0
-----

* Made minimal Blender requirement 2.80.
* Added `CHANGES.md` file.
* Pyrite now uses the N-panel (Sidebar) rather than the T-panel, as required
  for Blender 2.80 and above.
* Cleaned up the code (more PEP compliant).
* Start Over button also deletes Pyrite spheres.
* Minor improvements to user messages and plugin descriptions.

1.0.1
-----

* Trivial updates.

1.0.0
-----

* Original version.
