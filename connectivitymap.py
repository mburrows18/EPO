from numpy import sqrt, power

f1 = open('EPO_optimized.pdb', 'r')
f2 = open('epo.itp', 'r')
itpList = f2.read().splitlines()
pdbList = f1.read().splitlines()

#coordinates
coordinates = []
#connectivity
connectivity = []
#distance threshold
dist_thresh = 1.7

for line in pdbList:
    if line.startswith("HETATM"):
        words = line.split()                  # contains one line in the list format
        coordinates.append((words[1], (float(words[5]), float(words[6]), float(words[7]))))
#        print(coordinates)

for i in range(len(coordinates)):
    atom_coord_x, atom_coord_y, atom_coord_z = coordinates[i][1][0], coordinates[i][1][1], coordinates[i][1][2]
    connectivity.append([i + 1, []])
    for j in range(len(coordinates)):
        target_coord_x, target_coord_y, target_coord_z = coordinates[j][1][0], coordinates[j][1][1], coordinates[j][1][2]
        distance = sqrt(power((target_coord_x - atom_coord_x), 2) + 
            power((target_coord_y - atom_coord_y), 2) + 
            power((target_coord_z - atom_coord_z), 2))
        if (distance < dist_thresh) and (i != j):
            connectivity[i][1].append(j + 1)

print(connectivity)
            
            
