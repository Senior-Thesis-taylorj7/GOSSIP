import math as Math
import Bio.PDB
import EIGAs_utils
def characteristics(P1,q):
	ret = []
	C = centroid(P1)
	for a in range(len(P1)-q): #a is a residue
		temp1 = euclidean(P1[a], P1[a+2])
		temp2 = euclidean(P1[a+2], P1[a+q-2])
		temp3 = euclidean(P1[a+q-2],P1[a+q])
		temp4 = euclidean(P1[a+q], P1[a])
		temp5 = euclidean(P1[a], P1[a+q-2])
		P2 = centroid(P1)
		temp6 = Math.fabs(getAngle((P1[a][0] - P1[a+q][0], P1[a][1] - P1[a+q][1], P1[a][2] - P1[a+q][2]) , (P1[a][0] - P2[0], P1[a][1] - P2[1], P1[a][2] - P2[2]) ))
		temp7 = a
		temp8 = len(P1) - a + 1
		ret = ret + [[temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8]]
	return ret
def score(P1,P2):
	score = [[0 for x in range(len(P2))] for x in range(len(P1))] 
	for a in range(len(P1)):
		for b in range(len(P2)):
			score[a][b] = 3-2*Math.fabs(P1[a][3] - P2[b][3]) - 0.1*Math.fabs(P1[a][5] - P2[b][5])
	return score
def euclidean(a,b):
	return Math.pow(Math.pow(a[0]-b[0],2)+Math.pow(a[1]-b[1],2),0.5)
def centroid(list):
	x = 0.0
	y = 0.0
	z = 0.0
	for j in list:
		x += j[0]
		y += j[1]
		z += j[2]
	return ((x/len(list)),(y/len(list)),(z/len(list)))
def getAngle(a,b):
	return Math.acos((a[1]*b[1] + a[2]*b[2] + a[0]*b[0])/(Math.pow(a[1]*a[1]+a[0]*a[0]+a[2]*a[2],0.5)*Math.pow(b[1]*b[1]+b[0]*b[0] + b[2]*b[2],0.5)))
def normalize(a):
	max = None
	for item in a:
		size = Math.pow(item[0]*item[0] + item[1]*item[1] + item[2]*item[2],0.5)
		if(max == None or size > max):
			max = size
	ret = []
	for item in a:
		ret += [((item[0]+0.0)/max, (item[1]+0.0)/max,(item[2]+0.0)/max)]
	return ret
def itemDiv(sc):
	ret = []
	for i in sc:
		r2 = []
		for j in i:
			r2 += [j / max(sc)[0]]
		ret += [r2]
	return rets
logFile = open('../logs/runlog','w')
Prot = EIGAs_utils.readProteins(logFile, '../data/Skolnick40/', '../data/SCOP/dir.cla.scop.txt_1.75')
amk = Prot["1AMK"].model[0].chain[0].coords
aw2 = Prot["1AW2"].model[0].chain[0].coords
f = open("results.txt",'w')

#print score(characteristics(amk,15),characteristics(aw2,15))
for i in score(characteristics(amk,15),characteristics(aw2,15)):
	f.write(str(i))
f.close()

