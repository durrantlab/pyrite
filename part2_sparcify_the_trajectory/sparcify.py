import cStringIO as StringIO
import pymolecule
import numpy

# Here are some variables that eventually will be user specified.
sim_pdb_filename = "sim_to_test.pdb"
frame_stride = 2  # Only look at every other frame.
atom_stride = 5  # Only save every 5th atom.

# Use pymolecule to load in the pdb trajectory frames

# This list will store the frames as pymolecule objects
pdb_frames = []

current_frame = 1
first_line = True

# This list stores the atoms in the current frame
lines_in_this_frame = []

# Go through each of the lines in the multi-frame pdb file
for line in open(sim_pdb_filename):
	# Need to ignore first line?
	if first_line:
		first_line = False
	elif line.startswith("END"):
		# It's the end of the frame!
		if ((current_frame % frame_stride) == 1):
			# Make the pymolecule object.
			mol = pymolecule.Molecule()
			
			#Select only every certain number of atoms
			#print "From lines in frame: ", lines_in_this_frame[0:5]
			pdb_lines = numpy.array(lines_in_this_frame)
			#print "From all atoms: \n", pdb_lines[0:5]
			pdb_lines_atoms = pdb_lines[::atom_stride]
			#print "From selected atoms: \n", pdb_lines_atoms[0:5]

			# Load these pdb lines into that object
			final_pdb_lines = StringIO.StringIO("\n".join(pdb_lines_atoms))
			mol.load_pdb_into_using_file_object(final_pdb_lines, False, False, False)

			# Save this pymolecule object to our list
			pdb_frames.append(mol)
			print mol

		# Clear the list of lines for this frame, increment current frame
		lines_in_this_frame = []
		current_frame += 1
	else:
		# Add this line to the list of lines for this frame
		lines_in_this_frame.append(line)


# Ok, so pdb_frames is a list of pymolecule objects, each one corresponding to a frame.

for i, frame in enumerate(pdb_frames):
    	frame.save_pdb("frame" + str(i) + ".pdb", False, False, False)
