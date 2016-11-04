import cStringIO as StringIO
import pymolecule

# Here are some variables that eventually will be user specified.
sim_pdb_filename = "sim.pdb"
frame_stride = 2  # Only look at every other frame.
atom_stride = 5  # Only save every 5th atom.

# Use pymolecule to load in the pdb trajectory frames

# This list will store the frames as pymolecule objects
pdb_frames = []

# This list stores the atoms in the current frame
lines_in_this_frame = []

# Go through each of the lines in the multi-frame pdb file
for line in open(sim_pdb_filename):
	if line.startswith("END"):
		# It's the end of the frame!
		# Make the pymolecule object.
		mol = pymolecule.Molecule()

		# Load these pdb lines into that object
		pdb_lines = StringIO.StringIO("\n".join(lines_in_this_frame))
		mol.load_pdb_into_using_file_object(pdb_lines, False, False, False)

		# Save this pymolecule object to our list
		pdb_frames.append(mol)

		# Clear the list of lines for this frame
		lines_in_this_frame = []
		print mol

	# Add this line to the list of lines for this frame
	lines_in_this_frame.append(line)

# Ok, so pdb_frames is a list of pymolecule objects, each one corresponding to a frame.

for i, frame in enumerate(pdb_frames):
    	frame.save_pdb("frame" + str(i) + ".pdb", False, False, False)
