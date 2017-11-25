# Tested on Python 2.7.13 with MDAnalysis 0.16.2
# Examples of use:
#     python convert_binary_to_pdb.py my_coors.ncdf my_param.parm7
#     python convert_binary_to_pdb.py my_coors.dcd my_param.prmtop
#     python convert_binary_to_pdb.py my_pdb.pdb

from MDAnalysis import Universe, Writer
from MDAnalysis.analysis import align
import sys
import cStringIO
import os

# From the command line, the user will provide first the coordinate file, then
# the parameter file. The parameter file is option. If it is not specified,
# the coordinate file must include parameter information (e.g., the PDB
# format).
coord_file = sys.argv[1]
param_file = sys.argv[2] if len(sys.argv) > 2 else None

# Create the universe object. It contains the trajectory.
if param_file is None:
    u = Universe(coord_file)
else:
    u = Universe(param_file, coord_file)

# Get the first frame as PDB. Save it.
A = u.select_atoms("protein and name CA")
with Writer("__first_frame.tmp.pdb", A.n_atoms) as W:
    u.trajectory[0] # rewind trajectory
    W.write(A)

# Load it back in as the reference for a subsequent alignment.
ref = Universe("__first_frame.tmp.pdb")

# Align the trajectory to the reference (first) frame.
u.trajectory[0] # rewind trajectory
ref.trajectory[0]
alignment = align.AlignTraj(u, ref, select="protein and name CA",
                            weights='mass', in_memory=True)
alignment.run()

# Ignore atoms that belong to water molecules or counter ions. Retain only the
# hydrogen atoms on any ligands (not protein).
water_resnames = ["HOH", "WAT", "H2O", "TIP" ,"TIP3" ,"TP3", "TP4", "TIP3P"]
common_counter_ion_resnames = ["SOD", "POT", "CLA", "K+", "Na+", "Cl-"]
water_and_ions_resnames = water_resnames + common_counter_ion_resnames
water_and_ions_selection = 'resname ' + " ".join(water_and_ions_resnames)
not_water_ions_selection = '(not ' + water_and_ions_selection + ')'
no_Hs_selection = "(not name H*)"
lig_with_Hs_selection = "(" + not_water_ions_selection + " and (not protein))"
selection = "(" + not_water_ions_selection + " and " + no_Hs_selection + \
            ") or " + lig_with_Hs_selection
A = u.select_atoms(selection)

# To keep things simple, let's limit the max number of frames to 50. If there
# are more than 50 frames, start skipping frames.
max_frames = 50
skip = 1
while len(u.trajectory) / skip > max_frames:
    skip = skip + 1

# Save to the PDB file.
with Writer("output.pdb", A.n_atoms) as W:
    for ts in u.trajectory[::skip]:
        # Save the frame
        W.write(A)

# Remove the temp first-frame file
os.unlink("__first_frame.tmp.pdb")

# Let user know...
print("\nFile saved to output.pdb")