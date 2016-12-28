from numpy import sqrt, power, linalg
import numpy as np


f1 = open('EPO_optimized.pdb', 'r')
f2 = open('epo.itp', 'r')
f3 = open('connectivity.txt', 'w')
itpList = f2.read().splitlines()
pdbList = f1.read().splitlines()

#coordinates
coordinates = []
#connectivity by atom name
connectivity = []
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
#coord for angles in the form [xyz atom], [xyz central atom] [xyz atom]
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
    connectivitymap.append([i + 1, []])
    bdist.append([aname[i], [], []])
    atomcoords.append((aname[i], (atom_coord_x, atom_coord_y,atom_coord_z)))
    for j in range(len(coordinates)):
        target_coord_x, target_coord_y, target_coord_z = coordinates[j][1][0], coordinates[j][1][1], coordinates[j][1][2]
        distance = sqrt(power((target_coord_x - atom_coord_x), 2) + 
            power((target_coord_y - atom_coord_y), 2) + 
            power((target_coord_z - atom_coord_z), 2))
        if (distance < dist_thresh) and (i != j):
            connectivity[i][1].append(aname[j])
            connectivitymap[i][1].append(j + 1)
            bdist[i][1].append(j + 1)
            bdist[i][2].append(str(distance))
#print ("\n", file = f3)
print ("Bonding in EPO:", file = f3)
for i in range (len(connectivity)):
    #if len(connectivity[i][1]) == 1:
        #print ('{} is bonded to:'.format(aname[i]), file = f3)
        #print ('        {} \n with bond length {}'.format(', '.join(connectivity[i][1]), ', '.join(bdist[i][2])), file = f3)
        #print ('{} makes one bond'.format(aname[i]))
    if len(connectivity[i][1]) != 1:
        #print ('{} makes multiple bonds'.format(aname[i]))
        #print ('{} with coordinates {} is bonded to:'.format(aname[i], atomcoords[i][1]), file = f3)
        bondcoords.append([[aname[i], atomcoords[i][1]], []])
        bonds.append([aname[i], []])
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
            #print ('        {} with coordinates {}'.format(connectivity[i][1][j], atomcoords[Aindex][1]), file = f3)
            #print ('            {}-{} bond length = {} Angstroms'.format(aname[i], Aname, bdist[i][2][j]), file = f3)
            bondlengths.append((aname[i], Aname, bdist[i][2][j]))
            bondcoords[i][1].append([connectivity[i][1][j], atomcoords[Aindex][1]])
            #bondcoords[i][0] generates atom name with coordinates bondcoords[i][0][1]
            #bondcoords[i][1] generates [atom name, [coordinates]] of all bonded atoms
            bonds[i][1].append(connectivity[i][1][j])
   
for i in range (len(bonds)):
    print ('atom {} makes {} bonds \n    and is bonded to {}'.format(aname[i], len(connectivity[i][1]), ', '.join(bonds[i][1])), file = f3)   
    #print ('{} is bonded to {}'.format(bonds[i][0], ', '.join(bonds[i][1])), file = f3)
for i in range (len(bondcoords)):
    atomb = np.array(bondcoords[i][0][1])
    atoma = np.array(bondcoords[i][1][0][1])
    atomc = np.array(bondcoords[i][1][1][1])  
    F = atoma - atomb
    E = atomb - atomc
    dot = np.dot(F, E)
    normF = np.linalg.norm(F)
    normE = np.linalg.norm(E)
    cosangle = dot/(normE * normF)
    angle = 180 - np.degrees(np.arccos(cosangle))
    bondangle.append([[bondcoords[i][1][0][0], bondcoords[i][0][0], bondcoords[i][1][1][0]], []])
    anglenumber +=1
    numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][0][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][1][0])])])
    bondangle[anglenumber - 1][1].append(str(angle))
    if len(bondcoords[i][1]) == 3:
        atomb1 = np.array(bondcoords[i][0][1])
        atoma1 = np.array(bondcoords[i][1][1][1])
        atomc1 = np.array(bondcoords[i][1][2][1])
        F1 = atoma1 - atomb1
        E1 = atomb1 - atomc1
        dot1 = np.dot(F1, E1)
        normF1 = np.linalg.norm(F1)
        normE1 = np.linalg.norm(E1)
        cosangle1 = dot1/(normE1 * normF1)
        angle1 = 180 - np.degrees(np.arccos(cosangle1))
        bondangle.append([[bondcoords[i][1][1][0], bondcoords[i][0][0], bondcoords[i][1][2][0]], []])
        anglenumber +=1
        numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][1][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][2][0])])])
        bondangle[anglenumber - 1][1].append(str(angle1))
        atomb2 = np.array(bondcoords[i][0][1])
        atoma2 = np.array(bondcoords[i][1][0][1])
        atomc2 = np.array(bondcoords[i][1][2][1])
        F2 = atoma2 - atomb2
        E2 = atomb2 - atomc2
        dot2 = np.dot(F2, E2)
        normF2 = np.linalg.norm(F2)
        normE2 = np.linalg.norm(E2)
        cosangle2 = dot2/(normE2 * normF2)
        angle2 = 180 - np.degrees(np.arccos(cosangle2))
        bondangle.append([[bondcoords[i][1][0][0], bondcoords[i][0][0], bondcoords[i][1][2][0]], []])
        anglenumber +=1
        numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][0][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][2][0])])])
        bondangle[anglenumber - 1][1].append(str(angle2))
    if len(bondcoords[i][1]) == 4:
        atomb1 = np.array(bondcoords[i][0][1])
        atoma1 = np.array(bondcoords[i][1][1][1])
        atomc1 = np.array(bondcoords[i][1][2][1])
        F1 = atoma1 - atomb1
        E1 = atomb1 - atomc1
        dot1 = np.dot(F1, E1)
        normF1 = np.linalg.norm(F1)
        normE1 = np.linalg.norm(E1)
        cosangle1 = dot1/(normE1 * normF1)
        angle1 = 180 - np.degrees(np.arccos(cosangle1))
        bondangle.append([[bondcoords[i][1][1][0], bondcoords[i][0][0], bondcoords[i][1][2][0]], []])
        anglenumber +=1
        numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][1][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][2][0])])])
        bondangle[anglenumber - 1][1].append(str(angle1))
        atomb2 = np.array(bondcoords[i][0][1])
        atoma2 = np.array(bondcoords[i][1][0][1])
        atomc2 = np.array(bondcoords[i][1][2][1])
        F2 = atoma2 - atomb2
        E2 = atomb2 - atomc2
        dot2 = np.dot(F2, E2)
        normF2 = np.linalg.norm(F2)
        normE2 = np.linalg.norm(E2)
        cosangle2 = dot2/(normE2 * normF2)
        angle2 = 180 - np.degrees(np.arccos(cosangle2))
        bondangle.append([[bondcoords[i][1][0][0], bondcoords[i][0][0], bondcoords[i][1][2][0]], []])
        anglenumber +=1
        numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][0][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][2][0])])])
        bondangle[anglenumber - 1][1].append(str(angle2))
        atomb3 = np.array(bondcoords[i][0][1])
        atoma3 = np.array(bondcoords[i][1][0][1])
        atomc3 = np.array(bondcoords[i][1][3][1])
        F3 = atoma3 - atomb3
        E3 = atomb3 - atomc3
        dot3 = np.dot(F3, E3)
        normF3 = np.linalg.norm(F3)
        normE3 = np.linalg.norm(E3)
        cosangle3 = dot3/(normE3 * normF3)
        angle3 = 180 - np.degrees(np.arccos(cosangle3))
        bondangle.append([[bondcoords[i][1][0][0], bondcoords[i][0][0], bondcoords[i][1][3][0]], []])
        anglenumber +=1
        numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][0][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][3][0])])])
        bondangle[anglenumber - 1][1].append(str(angle3))
        atomb4 = np.array(bondcoords[i][0][1])
        atoma4 = np.array(bondcoords[i][1][2][1])
        atomc4 = np.array(bondcoords[i][1][3][1])
        F4 = atoma4 - atomb4
        E4 = atomb4 - atomc4
        dot4 = np.dot(F4, E4)
        normF4 = np.linalg.norm(F4)
        normE4 = np.linalg.norm(E4)
        cosangle4 = dot4/(normE4 * normF4)
        angle4 = 180 - np.degrees(np.arccos(cosangle4))
        bondangle.append([[bondcoords[i][1][2][0], bondcoords[i][0][0], bondcoords[i][1][3][0]], []])
        anglenumber +=1
        numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][2][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][3][0])])])
        bondangle[anglenumber - 1][1].append(str(angle4))
        atomb5 = np.array(bondcoords[i][0][1])
        atoma5 = np.array(bondcoords[i][1][1][1])
        atomc5 = np.array(bondcoords[i][1][3][1])
        F5 = atoma5 - atomb5
        E5 = atomb5 - atomc5
        dot5 = np.dot(F5, E5)
        normF5 = np.linalg.norm(F5)
        normE5 = np.linalg.norm(E5)
        cosangle5 = dot5/(normE5 * normF5)
        angle5 = 180 - np.degrees(np.arccos(cosangle5))
        bondangle.append([[bondcoords[i][1][1][0], bondcoords[i][0][0], bondcoords[i][1][3][0]], []])
        anglenumber +=1
        numberedangles.append([anglenumber, '-'.join([''.join(bondcoords[i][1][1][0]), ''.join(bondcoords[i][0][0]), ''.join(bondcoords[i][1][3][0])])])
        bondangle[anglenumber - 1][1].append(str(angle5))
for i in range (len(bondangle)):
    #print('angle number {}, {}, with bond angle = {} degrees'.format(numberedangles[i][0], '-'.join(bondangle[i][0]), ', '.join(bondangle[i][1])), file = f3)
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