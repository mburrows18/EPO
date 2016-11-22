f1 = open('EPO_optimized.pdb', 'r')
f2 = open('epo.itp', 'r')
itpList = f2.read().splitlines()
pdbList = f1.read().splitlines()
HETATMList = []
numberList = []
nameList = []
x_coordList = []
y_coordList = []
z_coordList = []
X_coordList = []
Y_coordList = []
Z_coordList = []
atypeList = []
massList = []
chargeList = []

for line in pdbList:
    if line.startswith("HETATM"):
       print ("yes") 
       list(line)
       HETATM = line[0:6]
       HETATMList.append(HETATM)
       HETATMList = [x.strip(' ') for x in HETATMList]
       number = line[9:11]
       numberList.append(number)
       numberList = [x.strip(' ') for x in numberList]
       name = line[12:16]
       nameList.append(name)
       nameList = [x.strip(' ') for x in nameList]
       x_coord = line[30:38]
       x_coordList.append(x_coord)
       x_coordList = [x.strip(' ') for x in x_coordList]
       X_coordList = [float(i) for i in x_coordList]
       y_coord = line[39:47]
       y_coordList.append(y_coord)
       y_coordList = [x.strip(' ') for x in y_coordList]
       Y_coordList = [float(i) for i in y_coordList]
       z_coord = line[48:55]
       z_coordList.append(z_coord)
       z_coordList = [x.strip(' ') for x in z_coordList]
       Z_coordList = [float(i) for i in z_coordList]
    else:
   	   print ("none")
print ("\n data from pdb file has been imported \n")

for line in itpList:
	if line.startswith(" "):
		print ("yes")
		list(line)
		atype = line[11:20]
		atypeList.append(atype)
		atypeList = [x.strip(' ') for x in atypeList]
		mass = line[59:68]
		massList.append(mass)
		massList = [x.strip(' ') for x in massList]
	else:
		print ("none")
print ("\n data from itp file has been imported \n")
print ("done")

import math
import numpy as np


from math import sqrt
for i in X_coordList:
	for j in Y_coordList:
		for k in Z_coordList:
			dx = X_coordList[0] - i
			dy = Y_coordList[0] - j
			dz = Z_coordList[0] - k
			distsq = pow(dx, 2) + pow(dy, 2) + pow(dz, 2)
			distance = sqrt(distsq)
			print (distance)

print ("final")
