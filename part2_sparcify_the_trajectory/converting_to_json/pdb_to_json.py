import json

num_frames = 50    # Set as number of the last frame you have stored as a pdb file in this folder
frame_list = []    # Will store lists of coordinate lists

# Go through each single-frame pdb file
for x in range(num_frames + 1):
    sim_pdb_filename = "frame" + str(x) + ".pdb"
    coordinate_list = []    # Will store lists of coordinates

    # Go through each of the lines in a single-frame pdb file
    for line in open(sim_pdb_filename):
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
    frame_list.append(coordinate_list)

# Write and save coordinate data in json file
with open('frame_coordinate_list.json', 'w') as f:
    json.dump(frame_list, f)