from numpy import sqrt, power
import numpy as np

f1 = open('EPO_optimized.pdb', 'r')
f2 = open('epo.itp', 'r')
f3 = open('connectivity.txt', 'w')
itpList = f2.read().splitlines()
pdbList = f1.read().splitlines()

#coordinates
coordinates = []
#connectivity by atom name without permutative duplicates
connectivity = []
#connectivity by atom name with duplicates
allconnectivity = []
#connectivity by atom index number
connectivitymap = []
#distance threshold
dist_thresh = 1.7
#atom type
atype = []
#atom name
aname = []
#bond lengths
bdist = []
#bonds in the form [atom, coordinates],  [[bonded atom, coordinates], [bonded atom, coordinates]], etc
bondcoords = []
#coord for angles in the form [xyz atom], [xyz central atom], [xyz atom]
anglecoords = []
#bonds in the form atom, [bonded atoms]
bonds = []
#all atoms with coords in the form [atom name, [(x, y, z)]]
atomcoords = []
#list of all bonds in form atom1, atom2, bond length
bondlengths = []
#bond angles in the form, [[bonded atom, central atom, bonded atom], bond angle]
bondangle = []
#angle number used to assign correct angle to bond
anglenumber = 0
#angle number with bond in the form [angle number, [atom1, central atom, atom2]]
numberedangles = []
#numerical id assigned to each angle
angleidnumbers = []
#angle name in the form atom1-central atom-atom2
angleidnames = []
#angles in degrees
anglevalues = []
index = -1

for line in itpList:
    if line.startswith(" "):
        words2 = line.split()
        atype.append(words2[1])
for line in pdbList:
    if line.startswith("HETATM"):
        words = line.split()                  # contains one line in the list format
        coordinates.append((words[1], (float(words[5]), float(words[6]), float(words[7]))))
        aname.append(words[2])
for i in range(len(coordinates)):
    atom_coord_x, atom_coord_y, atom_coord_z = coordinates[i][1][0], coordinates[i][1][1], coordinates[i][1][2]
    connectivity.append([aname[i], []])
    allconnectivity.append([aname[i], []])
    connectivitymap.append([i + 1, []])
    bdist.append([aname[i], [], []])
    atomcoords.append((aname[i], (atom_coord_x, atom_coord_y,atom_coord_z)))
    for j in range(len(coordinates)):
        target_coord_x, target_coord_y, target_coord_z = coordinates[j][1][0], coordinates[j][1][1], coordinates[j][1][2]
        distance = sqrt(power((target_coord_x - atom_coord_x), 2) + 
            power((target_coord_y - atom_coord_y), 2) + 
            power((target_coord_z - atom_coord_z), 2))
        if (distance < dist_thresh) and (i != j):
            allconnectivity[i][1].append(aname[j])
        if (distance < dist_thresh) and (i != j) and (i < j):
            connectivity[i][1].append(aname[j])
            connectivitymap[i][1].append(j + 1)
            bdist[i][1].append(j + 1)
            bdist[i][2].append(str(distance))
            bondcoords.append([[aname[i], atomcoords[i][1]], []])
            bonds.append([aname[i], []])
print ("Bonding in EPO:", file = f3)
for i in range (len(connectivity)):
    if len(connectivity[i][1]) !=0:
        print (connectivity[i], file = f3)
    if len(connectivity[i][1]) != 1:
        for j in range (len(connectivity[i][1])):
            #lines 76-84 are used to determine index of atom bonded to aname[i] in order to produce corresponding xyz coordinates
            Aname = connectivity[i][1][j]
            if Aname.startswith("C"):
                Aindex = (int(Aname.strip('C')) - 1)
            if Aname.startswith("H"):
                Aindex = (int(Aname.strip('H')) - 1)
            if Aname.startswith("N"):
                Aindex = (int(Aname.strip('N')) - 1)
            if Aname.startswith("O"):
                Aindex = (int(Aname.strip('O')) - 1)
            bondlengths.append((aname[i], Aname, bdist[i][2][j]))
            bondcoords[i][1].append([connectivity[i][1][j], atomcoords[Aindex][1]])
            #bondcoords[i][0] generates atom name with coordinates bondcoords[i][0][1]
            #bondcoords[i][1] generates [atom name, [coordinates]] of all bonded atoms
            bonds[i][1].append(connectivity[i][1][j])
for i in range (len(bonds)):
    if len(bondcoords[i][1]) !=0:
        atomb = np.array(bondcoords[i][0][1])
        atoma = np.array(bondcoords[i][1][0][1])
        atomc = np.array(bondcoords[i][1][1][1])  
        F = atoma - atomb
        E = atomb - atomc
        angle = np.dot(F, E)
        bondangle.append([[bondcoords[i][1][0][0], bondcoords[i][0][0], bondcoords[i][1][1][0]], []])
        anglenumber +=1
        numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][0][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][1][0])])])
        bondangle[anglenumber - 1][1].append(str(np.degrees(np.cos(angle))))
        if len(bondcoords[i][1]) == 3:
            atomb1 = np.array(bondcoords[i][0][1])
            atoma1 = np.array(bondcoords[i][1][1][1])
            atomc1 = np.array(bondcoords[i][1][2][1])
            F1 = atoma1 - atomb1
            E1 = atomb1 - atomc1
            angle1 = np.dot(F1, E1)
            bondangle.append([[bondcoords[i][1][1][0], bondcoords[i][0][0], bondcoords[i][1][2][0]], []])
            anglenumber +=1
            numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][1][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][2][0])])])
            bondangle[anglenumber - 1][1].append(str(np.degrees(np.cos(angle1))))
            atomb2 = np.array(bondcoords[i][0][1])
            atoma2 = np.array(bondcoords[i][1][0][1])
            atomc2 = np.array(bondcoords[i][1][2][1])
            F2 = atoma2 - atomb2
            E2 = atomb2 - atomc2
            angle2 = np.dot(F2, E2)
            bondangle.append([[bondcoords[i][1][0][0], bondcoords[i][0][0], bondcoords[i][1][2][0]], []])
            anglenumber +=1
            numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][0][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][2][0])])])
            bondangle[anglenumber - 1][1].append(str(np.degrees(np.cos(angle2))))
        if len(bondcoords[i][1]) == 4:
            atomb1 = np.array(bondcoords[i][0][1])
            atoma1 = np.array(bondcoords[i][1][1][1])
            atomc1 = np.array(bondcoords[i][1][2][1])
            F1 = atoma1 - atomb1
            E1 = atomb1 - atomc1
            angle1 = np.dot(F1, E1)
            bondangle.append([[bondcoords[i][1][1][0], bondcoords[i][0][0], bondcoords[i][1][2][0]], []])
            anglenumber +=1
            numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][1][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][2][0])])])                
            bondangle[anglenumber - 1][1].append(str(np.degrees(np.cos(angle1))))
            atomb2 = np.array(bondcoords[i][0][1])
            atoma2 = np.array(bondcoords[i][1][0][1])
            atomc2 = np.array(bondcoords[i][1][2][1])
            F2 = atoma2 - atomb2
            E2 = atomb2 - atomc2
            angle2 = np.dot(F2, E2)
            bondangle.append([[bondcoords[i][1][0][0], bondcoords[i][0][0], bondcoords[i][1][2][0]], []])
            anglenumber +=1
            numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][0][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][2][0])])])
            bondangle[anglenumber - 1][1].append(str(np.degrees(np.cos(angle2))))
            atomb3 = np.array(bondcoords[i][0][1])
            atoma3 = np.array(bondcoords[i][1][0][1])
            atomc3 = np.array(bondcoords[i][1][3][1])
            F3 = atoma3 - atomb3
            E3 = atomb3 - atomc3
            angle3 = np.dot(F3, E3)
            bondangle.append([[bondcoords[i][1][0][0], bondcoords[i][0][0], bondcoords[i][1][3][0]], []])
            anglenumber +=1
            numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][0][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][3][0])])])
            bondangle[anglenumber - 1][1].append(str(np.degrees(np.cos(angle3))))
            atomb4 = np.array(bondcoords[i][0][1])
            atoma4 = np.array(bondcoords[i][1][2][1])
            atomc4 = np.array(bondcoords[i][1][3][1])
            F4 = atoma4 - atomb4
            E4 = atomb4 - atomc4
            angle4 = np.dot(F4, E4)
            bondangle.append([[bondcoords[i][1][2][0], bondcoords[i][0][0], bondcoords[i][1][3][0]], []])
            anglenumber +=1
            numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][2][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][3][0])])])
            bondangle[anglenumber - 1][1].append(str(np.degrees(np.cos(angle4))))
            atomb5 = np.array(bondcoords[i][0][1])
            atoma5 = np.array(bondcoords[i][1][1][1])
            atomc5 = np.array(bondcoords[i][1][3][1])
            F5 = atoma5 - atomb5
            E5 = atomb5 - atomc5
            angle5 = np.dot(F5, E5)
            bondangle.append([[bondcoords[i][1][1][0], bondcoords[i][0][0], bondcoords[i][1][3][0]], []])
            anglenumber +=1
            numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][1][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][3][0])])])
            bondangle[anglenumber - 1][1].append(str(np.degrees(np.cos(angle5))))
for i in range (len(bondangle)):
    print ('# {}  {}  {}'.format(numberedangles[i][0], '-'.join(bondangle[i][0]), ''.join(bondangle[i][1])), file = f3)
f3.close()
f4 = open('connectivity.txt', 'r')
connectivityList = f4.read().splitlines()
for line in connectivityList:
    if line.startswith("#"):
        words3 = line.split()
        angleidnumbers.append(words3[1])
        angleidnames.append(words3[2])
        anglevalues.append(words3[3])      
f4.close()
f5 = open('connectivity.txt', 'w')
for line in connectivityList:
    if not line.startswith("#"):
        print(line, file = f5)
        index +=1
f5.close()
f5 = open('connectivity.txt', 'a')
for i in range (len(angleidnumbers)):
    print (angleidnumbers[i], angleidnames[i], anglevalues[i], file = f5)
print("All atom connectivity", file = f5)
for i in range (len(allconnectivity)):
    print (allconnectivity[i], file = f5)