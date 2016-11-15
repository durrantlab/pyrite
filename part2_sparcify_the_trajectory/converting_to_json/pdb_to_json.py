import json

num_frames = 20    # Set as number of the last frame you have stored as a pdb file in this folder
frame_list = []    # Will store lists of coordinate lists

# Go through each single-frame pdb file
for x in range(num_frames + 1):
    sim_pdb_filename = "frame" + str(num_frames) + ".pdb"
    coordinate_list = []    # Will store lists of coordinates

    # Go through each of the lines in a single-frame pdb file
    for line in open(sim_pdb_filename):
        # Obtain x, y, z coordinate elements from line if not empty
        if line != '\n':
            coordinates = [line[30] + line[31] + line[32] + line[33] + line[34] + line[35] + line[36] + line[37]]       # x
            coordinates.append(line[38] + line[39] + line[40] + line[41] + line[42] + line[43] + line[44] + line[45])   # y
            coordinates.append(line[46] + line[47] + line[48] + line[49] + line[50] + line[51] + line[52] + line[53])   # z
            coordinates = [float(coordinates[0].strip()), float(coordinates[1].strip()), float(coordinates[2].strip())]
            coordinate_list.append(coordinates)

    # Append list of coordinate lists from newly read frame to the frame list
    frame_list.append(coordinate_list)

# Write and save coordinate data in json file
with open('frame_coordinate_list.json', 'w') as f:
    json.dump(frame_list, f)