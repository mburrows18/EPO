
from numpy import sqrt, power, linalg
import numpy as np


f1 = open('EPO_optimized.pdb', 'r')
#input file
f2 = open('epo.itp', 'r')
#input file
f3 = open('connectivity.txt', 'w')
#intermediate output file

'''
atom number starts with 1
atom index starts with 0, so atom number 1 has atom index 0
'''

itpList = f2.read().splitlines()
#list of lines from itp input file, f2
pdbList = f1.read().splitlines()
#list of lines from pdb input file, f1
coordinates = []
#atom coordinates in the form [atom number, [x, y, z]]
atomcoords = []
#atom coordinates in the form [atom name, [x, y, z]]
connectivity = []
#bonded atoms by atom name in the form [atom1, [bonded atoms]]
uniqueconnectivity = []
#bonded atoms without permutative duplicates by atom name in the form [atom1, [bonded atoms]]
#in v1 of script this was also known as "bonds"
dist_thresh = 1.7
#distance threshold min bond length
aname = []
#atom names starting from C1 with index 0
atype = []
#atom type starting from atom type for C1 with index 0
bdist = []
#bond lengths in the form [atom1, [bonded atoms], [bond lenghts]]
uniquebdist = []
#bond lengths without permutative duplicates in the form [atom1, [bonded atoms], [bond lenghts]]
bondcoords = []
'''bond coordinates in the form 
		[[atom, coordinates], [[bonded atom, coordinates], etc]]
		bondcoords[i][0] generates atom name with coordinates bondcoords[i][0][1]
		bondcoords[i][1] generates [atom name, [coordinates]] of all bonded atoms'''
bondcoordtype = []
bondangletype = []
bondlengths = []
#all bond lengths with permutative duplicates in the form [atom1, atom2, bond length]
bondangle = []
#bond angles in the form [[bonded atom, central atom, bonded atom], bond angle]
numberedangles = []
#list of all angles in the form [angle number, [atom1, central atom, atom2]]
angleidnumbers = []
#stored list of angle numbers. angle number = index + 1 in this list
#extracted from connectivity.txt
angleidnames = []
#stored list of angles in the form [atom1-central atom-atom2] with angle number 1 having index 0
#extracted from connectivity.txt
angleidtypes = []
#stored list of angles in the form [atom1 type-central atom type-atom2 type] with angle number 1 having index 0
#extracted from connectivity.txt
angleIDtypes = []
anglevalues = []
#stored list of angle values with angle number 1 having index 0
#extracted from connectivity.txt
uniquebonds = []
#unique bonds in the form [atom1, atom2]. index = bond number - 1
uniquebondtypes = []
#unique bond types in the form [atomtype1, atomtype2] index = bond number - 1
uniquebondlengths = []
#unique bond lengths without permutative duplicates in the form [atom1, atom2, bond length]
bondtype = []
#unique bonds in the form [atom1 type, atom2 type]. index = index from unique bond list
Aindex = 0
#Aindex is NOT an ordered list -- intermediate variable used in section 2
anglenumber = 0
#anglenumber is NOT an ordered list -- intermediate variable used in section 3
uniquebondlengthsnm = []
#list of unique bond lengths in nanometers
KB = []
#list of bond KB values for bond types

def printlist(List):
#prints all items in list, output file = f3 (connectivity.txt)
	for i in range (len(List)):
		print (List[i], file = f3)
	return


def aindex(Aname):
#extracts index from atom name
	if Aname.startswith("C"):
		Aindex = (int(Aname.strip("C")) - 1)
	if Aname.startswith("H"):
		Aindex = (int(Aname.strip("H")) - 1)
	if Aname.startswith("N"):
		Aindex = (int(Aname.strip("N")) - 1)
	if Aname.startswith("O"):
		Aindex = (int(Aname.strip("O")) - 1)
	return Aindex


def calcangle(bondangle, numberedangles, anglenumber, bondcoords, bondcoordtype, i, a, c):
#calculates bond angles
	atomA = np.array(bondcoords[i][1][a][1])
	atomB = np.array(bondcoords[i][0][1])
	atomC = np.array(bondcoords[i][1][c][1])
	AB = atomA - atomB
	BC = atomB - atomC
	dot = np.dot(AB, BC)
	normAB = np.linalg.norm(AB)
	normBC = np.linalg.norm(BC)
	cosangle = dot/(normAB * normBC)
	angle = 180 - np.degrees(np.arccos(cosangle))
	anglenumber += 1
	bondangle.append([[
		bondcoords[i][1][a][0], 
		bondcoords[i][0][0], 
		bondcoords[i][1][c][0]], []])
	bondangletype.append([[
		bondcoordtype[i][1][a][0], 
		bondcoordtype[i][0][0], 
		bondcoordtype[i][1][c][0]]])
	numberedangles.append([anglenumber, '-'.join([
		''.join(bondcoords[i][1][a][0]), 
		''.join(bondcoords[i][0][0]), 
		''.join(bondcoords[i][1][c][0])])])
	bondangle[anglenumber - 1][1].append(str(angle))
	return anglenumber


'''
extracts data from: 
	itp file
	pdb file
generates: 
	coordinates[]
	connectivity[]
	uniqueconnectivity[]
	bdist[]
	uniquebdist[]
	uniquebonds[]
	atype[]
	aname[]
	atomcoords[]
'''
for line in itpList:
	if line.startswith(' '):
		words2 = line.split()
		atype.append(words2[1])
for line in pdbList:
	if line.startswith("HETATM"):
		words = line.split()
		coordinates.append((words[1], 
			(float(words[5]), float(words[6]), float(words[7]))))
		aname.append(words[2])
for i in range(len(coordinates)):
	atom_coord_x = coordinates[i][1][0]
	atom_coord_y = coordinates[i][1][1]
	atom_coord_z = coordinates[i][1][2]
	connectivity.append([aname[i], []])
	uniqueconnectivity.append([aname[i], []])
	bdist.append([aname[i], [], []])
	uniquebdist.append([aname[i], [], []])
	atomcoords.append((aname[i], 
		(atom_coord_x, atom_coord_y, atom_coord_z)))
	for j in range(len(coordinates)):
		target_coord_x = coordinates[j][1][0]
		target_coord_y = coordinates[j][1][1]
		target_coord_z = coordinates[j][1][2]
		distance = sqrt(
			power((target_coord_x - atom_coord_x), 2) +
			power((target_coord_y - atom_coord_y), 2) +
			power((target_coord_z - atom_coord_z), 2))
		if (distance < dist_thresh) and (i != j) and (i < j):
			uniquebonds.append([aname[i], aname[j]])
			uniquebondtypes.append([atype[i], atype[j]])
			uniqueconnectivity[i][1].append(aname[j])
			uniquebdist[i][1].append(aname[j])
			uniquebdist[i][2].append(str(distance))
		if (distance < dist_thresh) and (i != j):
			connectivity[i][1].append(aname[j])
			bdist[i][1].append(aname[j])
			bdist[i][2].append(str(distance))


'''
uses:
	aindex(Aname)
generates:
	bondcoords[]
	bondlengths[]
	uniquebondlengths[]
'''
for i in range(len(connectivity)):
	if len(connectivity[i][1]) != 0:
		bondcoords.append([[aname[i], atomcoords[i][1]], []])
		bondcoordtype.append([[atype[i], atomcoords[i][1]], []])
		for j in range(len(connectivity[i][1])):
			Aname = connectivity[i][1][j]
			Aindex = aindex(Aname)
			bondlengths.append((aname[i], Aname, bdist[i][2][j]))
			bondcoords[i][1].append([connectivity[i][1][j], atomcoords[Aindex][1]])
			bondcoordtype[i][1].append([atype[Aindex], atomcoords[Aindex][1]])
for i in range(len(uniqueconnectivity)):
	if len(uniqueconnectivity[i][1]) != 0:
		bondcoords.append([[aname[i], atomcoords[i][1]], []])
		for j in range(len(uniqueconnectivity[i][1])):
			Aname = uniqueconnectivity[i][1][j]
			Aindex = aindex(Aname)
			uniquebondlengths.append((aname[i], Aname, uniquebdist[i][2][j]))


'''
uses:
	calcangle(bondangle, numberedangles, anglenumber, bondcoords, i, a, c)
generates:
	bondangle[]
	numberedangles[]
'''
for i in range(len(bondcoords)):
	if len(bondcoords[i][1]) == 2:
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords,
			bondcoordtype, i, 0, 1)
	if len(bondcoords[i][1]) == 3:
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i, 0, 1)
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i, 1, 2)
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i, 0, 2)
	if len(bondcoords[i][1]) == 4:
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i, 0, 1)
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i, 1, 2)
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i,	0, 2)
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i, 0, 3)
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i, 2, 3)
		anglenumber = calcangle(
			bondangle, 
			numberedangles, 
			anglenumber, 
			bondcoords, 
			bondcoordtype, i, 1, 3) 
for i in range(len(bondangle)):
	print ('# {}	{}	{}'.format(
		numberedangles[i][0], 
		'-'.join(bondangle[i][0]), 
		''.join(bondangle[i][1])), 
		file = f3)
	print ('2# {}	{}'.format(
		numberedangles[i][0], 
		'-'.join(bondangletype[i][0])), 
		file = f3)
print ("unique bonds", file = f3)
for i in range (len(uniquebonds)):
	print (i+1,uniquebonds[i], file = f3)
print ("unique bond lengths", file = f3)
for i in range(len(uniquebondlengths)):
	print (i+1,uniquebondlengths[i], file = f3)
print ("atom types", file = f3)
for i in range(len(aname)):
	print ('{} -- {}'.format(aname[i], atype[i]), file = f3)
f3.close()

f4 = open('connectivity.txt', 'r')
connectivityList = f4.read().splitlines()
for line in connectivityList:
  if line.startswith("#"):
    words3 = line.split()
    angleidnumbers.append(words3[1])
    angleidnames.append(words3[2])
    anglevalues.append(words3[3])
  if line.startswith("2#"):
    words4 = line.split()
    angleidtypes.append(words4[2])    
f4.close()
f5 = open('connectivity.txt', 'w')
for line in connectivityList:
  if not line.startswith("#") and not line.startswith("2#"):
    print(line, file = f5)
f5.close()
f5 = open('connectivity.txt', 'a')
print ("angles", file = f5)
for i in range (len(angleidnumbers)):
	print (angleidnumbers[i], 
		angleidnames[i], 
		angleidtypes[i], 
		anglevalues[i], 
		file = f5)
f5.close()


f6 = open('EPO.prm', 'w')
print ("\n[ bondtypes ]", file = f6)
print(";", 
	str("i").ljust(8),
	str("j").ljust(8),
	str("func").ljust(5),
	str("b0").ljust(25),
	str("kb").ljust(12),
	file = f6)
for i in range(len(uniquebondlengths)):
	lengthnm = float(uniquebondlengths[i][2])*0.1
	uniquebondlengthsnm.append([lengthnm])
typei = 0
typej = 0
kb = []

for i in range(len(uniquebondtypes)):
	if uniquebondtypes[i][0] == ("CN3") and uniquebondtypes[i][1] == ("CN3"):
		KB.append("418400.00")
	if uniquebondtypes[i][0] == ("CN3") and uniquebondtypes[i][1] == ("HN3"):
		KB.append("292880.00")
	if uniquebondtypes[i][0] == ("CN2") and uniquebondtypes[i][1] == ("CN3"):
		KB.append("267776.00")
	if uniquebondtypes[i][0] == ("CN3") and uniquebondtypes[i][1] == ("CN2"):
		KB.append("267776.00")
	if uniquebondtypes[i][0] == ("CN2") and uniquebondtypes[i][1] == ("CN2"):
		KB.append("252713.60")
	if uniquebondtypes[i][0] == ("CN2") and uniquebondtypes[i][1] == ("OG3R60"):
		KB.append("UNKNOWN")
	if uniquebondtypes[i][0] == ("CN2") and uniquebondtypes[i][1] == ("NN2"):
		KB.append("334720.00")
	if uniquebondtypes[i][0] == ("NN2") and uniquebondtypes[i][1] == ("CN2"):
		KB.append("334720.00")
	if uniquebondtypes[i][0] == ("NN2") and uniquebondtypes[i][1] == ("CN8"):
		KB.append("334720.00")
	if uniquebondtypes[i][0] == ("OG3R60") and uniquebondtypes[i][1] == ("CN2"):
		KB.append("UNKNOWN")
	if uniquebondtypes[i][0] == ("CN8") and uniquebondtypes[i][1] == ("CN7"):
		KB.append("186188.00")
	if uniquebondtypes[i][0] == ("CN8") and uniquebondtypes[i][1] == ("HN8"):
		KB.append("258571.20")
	if uniquebondtypes[i][0] == ("CN7") and uniquebondtypes[i][1] == ("HN7"):
		KB.append("258571.20")

for i in range(len(uniquebonds)):
	print(" ", 
		str(uniquebondtypes[i][0]).ljust(8), 
		str(uniquebondtypes[i][1]).ljust(8), 
		str("1").ljust(5), 
		str(uniquebondlengthsnm[i][0]).ljust(25), 
		str(KB[i]).ljust(12), 
		file = f6)
print ("\n\n", file = f6)
print ("[ angletypes ]", file = f6)
print(";", 
	str("i").ljust(8),
	str("j").ljust(8),
	str("k").ljust(8),
	str("func").ljust(5),
	str("theta0").ljust(15),
	str("ktheta").ljust(8),
	str("ub0").ljust(8),
	str("kub").ljust(8), 
	file = f6)
for i in range(len(angleidtypes)):
	line = str(angleidtypes[i])
	angleIDtypes.append(line.split("-"))

for i in range(len(angleidtypes)):
	print(" ", 
		str(angleIDtypes[i][0]).ljust(8), 
		str(angleIDtypes[i][1]).ljust(8), 
		str(angleIDtypes[i][2]).ljust(8), 
		str("5").ljust(5),
		str(anglevalues[i]).ljust(15), 
		#str(KTHETA).ljust(8),
		#str(UB0).ljust(8),
		#str(KUB).ljust(8), 
		file = f6)




