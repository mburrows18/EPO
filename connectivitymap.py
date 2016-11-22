f1 = open('EPO_optimized.pdb', 'r')
f2 = open('epo.itp', 'r')
itpList = f2.read().splitlines()
pdbList = f1.read().splitlines()
HETATMList = []
#atom number
numberList = []
#atom name
nameList = []
#x coordinates
x_coordList = []
X_coordList = []
#y coordinates
y_coordList = []
Y_coordlist = []
#z coordinates
z_coordList = []
Z_coordList = []
#atom type
atypeList = []
#atom mass
massList = []
#atom charge
chargeList = []
#interatomic distance
adistance = []
Adistance = []
#interatomic distances < 1.7
bdistance = []
#zip list of interacting atoms
Pairs = []
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


'''dx = (X_coordList[0] - X_coordList[1])
dy = (Y_coordList[0] - Y_coordList[1])
dz = (Z_coordList[0] - Z_coordList[1])
distsq = pow(dx, 2) + pow(dy, 2) + pow(dz, 2)
distance = sqrt(distsq)
print (distance)

print ("Test Distance")
'''

for a in range (0, 29, 1):
  Xcoord = X_coordList[a]
  Ycoord = Y_coordList[a]
  Zcoord = Z_coordList[a]
  #print (a)
  a1 = a
  for i in range (0, 29, 1):
    xcoord = X_coordList[i]
    ycoord = Y_coordList[i]
    zcoord = Z_coordList[i]
    dx = (xcoord - Xcoord)
    dy = (ycoord - Ycoord)
    dz = (zcoord - Zcoord)
    distsq = pow(dx, 2) + pow(dy, 2) + pow(dz, 2)
    distance = sqrt(distsq)
    adistance.append(distance)
    a2 = i
    print (a1, a2)
        #print (i, distance)
Adistance = [float(z) for z in adistance]
for d in Adistance:
  if d < 1.5:
    bdistance.append(d)
  else:
    print ("nonbonded")
print (bdistance)
# bdistance contains lengths of all "bonds"
print ("final")
