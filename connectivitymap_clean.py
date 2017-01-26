
from numpy import sqrt, power, linalg, average
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
uniqueatypes = []
#list of unique atom types in molecule
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
#same as bondcoords[] with atom names replaced with atom types
bondlengths = []
#all bond lengths with permutative duplicates in the form [atom1, atom2, bond length]
bondangle = []
#bond angles in the form [[bonded atom, central atom, bonded atom], bond angle]
bondangletype = []
#same as bondangle[] with atom names replaced with atom types
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
#stored list of angles in the form [[atom1 type], [central atom type], [atom2 type]] with angle number 1 having index 0
#xtracted from connectivity.txt
anglevalues = []
#stored list of angle values with angle number 1 having index 0
#extracted from connectivity.txt
uniquebonds = []
#unique bonds in the form [atom1, atom2]. index = bond number - 1
uniquebondtypes = []
#unique bond types in the form [atomtype1, atomtype2] index = bond number - 1
uniquebondlengths = []
#unique bond lengths without permutative duplicates in the form [atom1, atom2, bond length]
Aindex = 0
#Aindex is NOT an ordered list -- intermediate variable used in section 2
anglenumber = 0
#anglenumber is NOT an ordered list -- intermediate variable used in section 3
uniquebondlengthsnm = []
#list of unique bond lengths in nanometers
KB = []
#list of bond KB values for bond types
KBA = []
#list of angle KB values for angle types
bondtypelist = []
#stores data extracted from raw .prm file
uniquebtypelist = []
#intermediate list in removing permutative bonds
uniquebondtypelist = []
#final list of unique bond types and averaged bond lengths, prints to final .prm file
uniquebtype = []
#intermediate list in removing permutative bonds
btypelist = []
#intermediate list in removing permutative bonds in the form bondtypenum, atom1typenum, atom2typenum
avgblengths = []
#intermediate list in removing permutative bonds
blengths = []
#intermediate list in removing permutative bonds

#prints all items in list, output file = f3 (connectivity.txt)
def printlist(List):
	for i in range (len(List)):
		print (List[i])
	return

#averages all values in list LIST, returns avg as average
def average(LIST):
	avg = np.average(LIST)
	return avg

#extracts index from atom name
def aindex(Aname):
	if Aname.startswith("C"):
		Aindex = (int(Aname.strip("C")) - 1)
	if Aname.startswith("H"):
		Aindex = (int(Aname.strip("H")) - 1)
	if Aname.startswith("N"):
		Aindex = (int(Aname.strip("N")) - 1)
	if Aname.startswith("O"):
		Aindex = (int(Aname.strip("O")) - 1)
	return Aindex


#calculates bond angles
def calcangle(bondangle, numberedangles, anglenumber, bondcoords, bondcoordtype, i, a, c):
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
	uniquebondtypes[]
	atype[]
	uniqueatypes[]
	aname[]
	atomcoords[]
	words2
	words
'''
for line in itpList:
	if line.startswith(' '):
		words2 = line.split()
		atype.append(words2[1])
for i in atype:
	if i not in uniqueatypes:
		uniqueatypes.append(i)		
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
	uniqueconnectivity[]
	Aname
	Aindex
	connectivity[]
generates:
	bondcoords[]
	bondcoordtype[]
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
	uniquebonds[]
	uniquebondlengths[]
generates:
	bondangle[]
	numberedangles[]
  bondangletype[]
prints data to:
	connectivity.txt (f3)
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


'''
extracts data from:
	connectivity.txt (f4)
generates:
	angleidnumbers[]
	angleidnames[]
	angleidtypes[]
	anglevalues[]
	connectivityList
	words3
	words4
prints data to:
	connectivity.txt (f5)
'''
f4 = open('connectivity.txt', 'r')
#intermediate output file
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
#intermediate output file
for line in connectivityList:
  if not line.startswith("#") and not line.startswith("2#"):
    print(line, file = f5)
f5.close()
f5 = open('connectivity.txt', 'a')
#intermediate output file
print ("angles", file = f5)
for i in range (len(angleidnumbers)):
	print (angleidnumbers[i], 
		angleidnames[i], 
		angleidtypes[i], 
		anglevalues[i], 
		file = f5)
f5.close()


'''
uses:
	uniquebondlengths[]
	uniquebondtypes[]
	uniquebonds[]
	angleidtypes[]
	anglevalues[]
generates:
	uniquebondlengthsnm[]
	angleIDtypes[]
	KB[]
	KBA[]
prints data to:
	EPO.prm (f6)
'''
f6 = open('EPO.prm', 'w')
#.prm file for GROMACS simulations
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
	print("  ", 
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
	str("ktheta").ljust(12),
	str("ub0").ljust(10),
	str("kub").ljust(8), 
	file = f6)
for i in range(len(angleidtypes)):
	line = str(angleidtypes[i])
	angleIDtypes.append(line.split("-"))
for i in range(len(angleIDtypes)):
	if angleIDtypes[i][0] == ("CN3") and angleIDtypes[i][1] == ("CN3") and angleIDtypes[i][2] == ("CN3"):
		KBA.append("334.720000")
	if angleIDtypes[i][0] == ("CN3") and angleIDtypes[i][1] == ("CN3") and angleIDtypes[i][2] == ("HN3"):
		KBA.append("317.984000")
	if angleIDtypes[i][0] == ("CN3") and angleIDtypes[i][1] == ("CN3") and angleIDtypes[i][2] == ("CN2"):
		KBA.append("711.280000")
	if angleIDtypes[i][0] == ("CN3") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("CN2"):
		KBA.append("502.080000")
	if angleIDtypes[i][0] == ("CN2") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("OG3R60"):
		KBA.append("UNKNOWN")
	if angleIDtypes[i][0] == ("CN3") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("OG3R60"):
		KBA.append("UNKNOWN")
	if angleIDtypes[i][0] == ("CN2") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("CN3"):
		KBA.append("502.080000")
	if angleIDtypes[i][0] == ("CN3") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("NN2"):
		KBA.append("836.800000")
	if angleIDtypes[i][0] == ("CN2") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("NN2"):
		KBA.append("836.800000")
	if angleIDtypes[i][0] == ("CN2") and angleIDtypes[i][1] == ("CN3") and angleIDtypes[i][2] == ("HN3"):
		KBA.append("317.984000")
	if angleIDtypes[i][0] == ("CN2") and angleIDtypes[i][1] == ("CN3") and angleIDtypes[i][2] == ("CN3"):
		KBA.append("711.280000")
	if angleIDtypes[i][0] == ("CN2") and angleIDtypes[i][1] == ("NN2") and angleIDtypes[i][2] == ("CN2"):
		KBA.append("753.120000")
	if angleIDtypes[i][0] == ("CN2") and angleIDtypes[i][1] == ("NN2") and angleIDtypes[i][2] == ("CN8"):
		KBA.append("585.760000")
	if angleIDtypes[i][0] == ("CN2") and angleIDtypes[i][1] == ("OG3R60") and angleIDtypes[i][2] == ("CN2"):
		KBA.append("UNKNOWN")
	if angleIDtypes[i][0] == ("NN2") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("CN2"):
		KBA.append("836.800000")
	if angleIDtypes[i][0] == ("NN2") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("CN3"):
		KBA.append("836.800000")
	if angleIDtypes[i][0] == ("OG3R60") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("CN2"):
		KBA.append("UNKNOWN")
	if angleIDtypes[i][0] == ("OG3R60") and angleIDtypes[i][1] == ("CN2") and angleIDtypes[i][2] == ("CN3"):
		KBA.append("UNKNOWN")
	if angleIDtypes[i][0] == ("NN2") and angleIDtypes[i][1] == ("CN8") and angleIDtypes[i][2] == ("CN7"):
		KBA.append("920.480000")
	if angleIDtypes[i][0] == ("CN7") and angleIDtypes[i][1] == ("CN8") and angleIDtypes[i][2] == ("HN8"):
		KBA.append("288.947040")
	if angleIDtypes[i][0] == ("NN2") and angleIDtypes[i][1] == ("CN8") and angleIDtypes[i][2] == ("CN8"):
		KBA.append("UNKNOWN")
	if angleIDtypes[i][0] == ("HN8") and angleIDtypes[i][1] == ("CN8") and angleIDtypes[i][2] == ("HN8"):
		KBA.append("288.947040")
	if angleIDtypes[i][0] == ("CN8") and angleIDtypes[i][1] == ("CN7") and angleIDtypes[i][2] == ("HN7"):
		KBA.append("288.696000")
	if angleIDtypes[i][0] == ("HN7") and angleIDtypes[i][1] == ("CN7") and angleIDtypes[i][2] == ("HN7"):
		KBA.append("334.720000")
	if angleIDtypes[i][0] == ("NN2") and angleIDtypes[i][1] == ("CN8") and angleIDtypes[i][2] == ("HN8"):
		KBA.append("279.742240")
for i in range(len(angleidtypes)):
	print(" _", 
		str(angleIDtypes[i][0]).ljust(8), 
		str(angleIDtypes[i][1]).ljust(8), 
		str(angleIDtypes[i][2]).ljust(8), 
		str("5").ljust(5),
		str(anglevalues[i]).ljust(15), 
		str(KBA[i]).ljust(8),
		str(" -- ").ljust(8),
		str(" -- ").ljust(8), 
		file = f6)
f6.close()



#FROM HERE DOWN THE CODE WORKS, BUT IS QUITE MESSY. ALL CODE ABOVE THIS LINE HAS BEEN OPTIMIZED AND SIMPLIFIED AS MUCH AS POSSIBLE


'''
extracts data from:
	EPO.prm (f7)
generates:
	bondtypelist = []
	angletypelist = []
	uniquebtypelist = []
	uniquebondtypelist = []
	uniquebtype = []
	btypelist = []
	avgblengths = []
	blengths = []
'''
angletypelist = []
uniqueatypelist = []
atypelist = []
uniqueangletypelist = []
avgatheta = []
atheta = []
uniqueatype = []

f7 = open('EPO.prm', 'r')
prmList = f7.read().splitlines()
for line in prmList:
	if line.startswith(" _ "):
		linea = line.split(" _ ", 1)[1]
		angletypeList = linea.split()
		for i in range(len(uniqueatypes)):
			if angletypeList[0] == uniqueatypes[i]:
				Type1num = i
			if angletypeList[1] == uniqueatypes[i]:
				Type2num = i
			if angletypeList[2] == uniqueatypes[i]:
				Type3num = i
		angletypelist.append(
			[[], str(angletypeList[0]), Type1num, 
			str(angletypeList[1]), Type2num, 
			str(angletypeList[2]), Type3num, 
			str(angletypeList[3]), 
			str(angletypeList[4]), 
			str(angletypeList[5]), 
			str(angletypeList[6]), 
			str(angletypeList[7]), []])
	if line.startswith("  "):
		bondtypeList = line.split()
		for i in range(len(uniqueatypes)):
			if bondtypeList[0] == uniqueatypes[i]:
				Type1num = i
			if bondtypeList[1] == uniqueatypes[i]:
				Type2num = i
		bondtypelist.append(
			[[], str(bondtypeList[0]), Type1num, 
			str(bondtypeList[1]), Type2num, 
			str(bondtypeList[2]), 
			str(bondtypeList[3]), 
			str(bondtypeList[4]), []])
for i in range(len(angletypelist)):
	atom1T = angletypelist[i][2]
	atom2T = angletypelist[i][4]
	atom3T = angletypelist[i][6]
	angleT = ('{}-{}-{}'.format(atom1T, atom2T, atom3T))
	if angleT not in uniqueatypelist:
		uniqueatypelist.append(angleT)
		atypelist.append([[], atom1T, atom2T, atom3T])
for i in range(len(atypelist)):
	atypelist[i][0] = str(i)
	for j in range(len(atypelist)):
		atheta.append([j, []])
		avgatheta.append([])
		if atypelist[i][1] == atypelist[j][3] and atypelist[i][3] == atypelist[j][1] and atypelist[i][2] == atypelist[j][2] and (i != j):
			atypelist[i][0] = str('- {} {}'.format(j, i))
			atypelist[j][0] = str('- {} {}'.format(j, i))
for j in range(len(atypelist)):
	for i in range(len(angletypelist)):
		if not atypelist[j][0].startswith("-") and angletypelist[i][2] == atypelist[j][1] and angletypelist[i][4] == atypelist[j][2] and angletypelist[i][6] == atypelist[j][3]:
			atheta[j][1].append(float(angletypelist[i][8]))
			angletypelist[i][0] = str(j)
		if atypelist[j][0].startswith("-") and angletypelist[i][2] == atypelist[j][1] and angletypelist[i][4] == atypelist[j][2] and angletypelist[i][6] == atypelist[j][3] or atypelist[j][0].startswith("-") and angletypelist[i][2] == atypelist[j][3] and angletypelist[i][4] == atypelist[j][2] and angletypelist[i][6] == atypelist[j][1]:
			atheta[j][1].append(float(angletypelist[i][8]))
			angletypelist[i][0] = str(j)
	avgatheta[j] = average(atheta[j][1])
for i in range(len(angletypelist)):
	for j in range(len(avgatheta)):
		if angletypelist[i][0] == str(j) and avgatheta[j] not in angletypelist[i][12]:
			angletypelist[i][12] = avgatheta[j]
for i in range(len(angletypelist)):
	angletype = angletypelist[i][0]
	if angletype not in uniqueatype:
		uniqueatype.append(angletype)
		uniqueangletypelist.append(
			[str(angletypelist[i][0]), 
			str(angletypelist[i][1]), 
			str(angletypelist[i][2]), 
			str(angletypelist[i][3]),
			str(angletypelist[i][4]),
			str(angletypelist[i][5]),
			str(angletypelist[i][6]),
			str(angletypelist[i][7]),
			str(angletypelist[i][8]),
			str(angletypelist[i][9]),
			str(angletypelist[i][10]),
			str(angletypelist[i][11]),
			str(angletypelist[i][12])])
for i in range(len(bondtypelist)):
	atom1T = bondtypelist[i][2]
	atom2T = bondtypelist[i][4]
	bondT = ('{}-{}'.format(atom1T, atom2T))
	if bondT not in uniquebtypelist:
		uniquebtypelist.append(bondT)
		btypelist.append([[], atom1T, atom2T])
for i in range(len(btypelist)):
	btypelist[i][0] = str(i)
	for j in range(len(btypelist)):
		blengths.append([j, []])
		avgblengths.append([])
		if btypelist[i][1] == btypelist[j][2] and btypelist[i][2] == btypelist[j][1] and (i != j):
			btypelist[i][0] = str('- {} {}'.format(j, i))
			btypelist[j][0] = str('- {} {}'.format(j, i))
for j in range(len(btypelist)):
	for i in range(len(bondtypelist)):
		if not btypelist[j][0].startswith("-") and bondtypelist[i][2] == btypelist[j][1] and bondtypelist[i][4] == btypelist[j][2]:
			blengths[j][1].append(float(bondtypelist[i][6]))
			bondtypelist[i][0] = str(j)
		if btypelist[j][0].startswith("-") and bondtypelist[i][2] == btypelist[j][1] and bondtypelist[i][4] == btypelist[j][2] or btypelist[j][0].startswith("-") and bondtypelist[i][2] == btypelist[j][2] and bondtypelist[i][4] == btypelist[j][1]:
			blengths[j][1].append(float(bondtypelist[i][6]))
			bondtypelist[i][0] = str(j)
	avgblengths[j] = average(blengths[j][1])
for i in range(len(bondtypelist)):
	for j in range(len(avgblengths)):	
		if bondtypelist[i][0] == str(j) and avgblengths[j] not in bondtypelist[i][8]:
			bondtypelist[i][8] = avgblengths[j]
for i in range(len(bondtypelist)):
	bondtype = bondtypelist[i][0]
	if bondtype not in uniquebtype:
		uniquebtype.append(bondtype)
		uniquebondtypelist.append(
			[str(bondtypelist[i][0]), 
			str(bondtypelist[i][1]), 
			str(bondtypelist[i][2]),
			str(bondtypelist[i][3]), 
			str(bondtypelist[i][4]), 
			str(bondtypelist[i][5]),
			str(bondtypelist[i][6]), 
			str(bondtypelist[i][7]), 
			str(bondtypelist[i][8])])
f7.close()


'''
uses:
	uniquebondtypelist[]
	prmList[]
prints data to:
	EPO.prm (f8)
'''
f8 = open('EPO.prm', 'w')
for i in range(0, 3):
	print(prmList[i], file = f8)
for i in range(len(uniquebondtypelist)):
	print(" ", 
		str(uniquebondtypelist[i][1]).ljust(8), 
		str(uniquebondtypelist[i][3]).ljust(8), 
		str(uniquebondtypelist[i][5]).ljust(5), 
		str(uniquebondtypelist[i][8]).ljust(23), 
		str(uniquebondtypelist[i][7]).ljust(12), 
		file = f8)
for i in range(34, 39):
	print(prmList[i], file = f8)
for i in range(len(uniqueangletypelist)):
	print(" ",
		str(uniqueangletypelist[i][1]).ljust(8), 
		str(uniqueangletypelist[i][3]).ljust(8), 
		str(uniqueangletypelist[i][5]).ljust(8), 
		str(uniqueangletypelist[i][7]).ljust(5),
		str(uniqueangletypelist[i][12]).ljust(15), 
		str(uniqueangletypelist[i][9]).ljust(12),
		str(uniqueangletypelist[i][10]).ljust(10),
		str(uniqueangletypelist[i][11]).ljust(8), 
		file = f8)