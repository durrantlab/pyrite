import cStringIO as StringIO
import scoria
import numpy
import json

class Mineral:
    def __init__(self):
        self.pdb_filename = ""
        self.frame_stride = 1
        self.atom_stride = 1
        self.sparce_pdb_txt_frames = []
        self.json_rep = ""
        

    def load_pdb_trajectory(self, pdb_filename, frame_stride, atom_stride):
        self.pdb_filename = pdb_filename
        self.frame_stride = frame_stride
        self.atom_stride = atom_stride

        # # Here are some variables that eventually will be user specified.
        # sim_pdb_filename = "M2_traj.pdb"
        # frame_stride = 2  # Only look at every other frame.
        # atom_stride = 50  # Only save every 5th atom.

        # Load the trajectory
        traj = scoria.Molecule()
        traj.load_pdb_trajectory_into(self.pdb_filename, bonds_by_distance = False, serial_reindex = False, resseq_reindex = False)

        # Select every few atoms for now. But when we have the code to make certain regions denser, this will get more complicated.
        sel = traj.select_all()[::atom_stride]

        # Make a new molecule object from that selection, given the current frame.
        pruned_traj = traj.get_molecule_from_selection(sel, serial_reindex = False, resseq_reindex = False)

        # Save only every few frames (the pdb txt)
        for frame_index in range(0, traj.get_trajectory_frame_count(), self.frame_stride):
            txt = pruned_traj.save_pdb("", False, False, True, frame_index)
            self.sparce_pdb_txt_frames.append(txt)
    
    def make_json_representation(self):
        frames_list = []    # Will store lists of coordinate lists

        # Go through each single-frame pdb file
        for x in range(len(self.sparce_pdb_txt_frames)):
            coordinate_list = []    # Will store lists of coordinates

            txt = self.sparce_pdb_txt_frames[x].strip().split("\n")
            # Go through each of the lines in a single-frame pdb file
            for line in txt:
                # Obtain x, y, z coordinate elements from line if not empty
                if line != '\n':
                    coordinates = [line[30:38]]       # x
                    coordinates.append(line[38:46])   # y
                    coordinates.append(line[46:54])   # z
                    coordinates = [float(coordinates[0].strip()), float(coordinates[1].strip()), float(coordinates[2].strip())]

                    coordinate_list.append(coordinates)

                    # 0, 1, 2
                    ## 0, 2, 1
                    ## 1, 0, 2
            # Append list of coordinate lists from newly read frame to the frame list
            frames_list.append(coordinate_list)

        # Write and save coordinate data in json file
        # In the future, everything will be in one big class:
        # self.json_rep =json.dumps(frames_list)

        # For now, you need a way of loading into blender, so make a separate json file.
        with open('frame_coordinate_list.json', 'w') as f:
            json.dump(frames_list, f)        

    
mineral = Mineral()
mineral.load_pdb_trajectory("M2_traj.pdb", 2, 50)
mineral.make_json_representation()
