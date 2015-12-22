import math
import numpy as np

th_sim = 1
th_dist = 10

def dist(c1,c2):
    return math.pow(sum([(c1[i]-c2[i])**2 for i in range(len(c1))]),.5)

def Get_Angle(c0,c1,c2):
    P12 = dist(c0,c1)
    P13 = dist(c0,c2)
    P23 = dist(c1,c2)
    return np.arccos((P12**2 + P13**2 - P23**2) / (2 * P12 * P13))

def Get_Quadrilateral_Signatures(C,quadsize):
    qs = []
    centroid = [0.0,0.0,0.0]
    for i in range(0,len(C)):
        centroid = np.add(centroid, C[i])
    centroid = np.divide(centroid,len(C))
    for i in range(0,len(C)-quadsize):
        char = []
        coords = [C[i], C[i+2],C[i+quadsize-2],C[i+quadsize]]
        for j in range(0, len(coords)):
            char += [dist(coords[j], coords[j%4])]#Side Lengths
        char += [dist(coords[0], coords[2])]#Diagonals
        char += [get_angle(coords[0], coords[3], centroid)]#Theta
        char += [i+1,len(C)-i] #no +1 because 0 indexing
        qs += [char]
    return qs

Side_Resolution = 1
Diag_Resolution = 2*Side_Resolution
small_Angle = math.pi*25.0/180.0
large_Angle = math.pi*15.0/180.0
medium_Angle = math.pi*20.0/180.0

def get_Ind_Resolution(Prot1,Prot2):
    return (1.0-th_sim)*max(len(Prot1),len(Prot2))

def get_Quad_Size(Prot1,Prot2):
    size = max(len(Prot1),len(Prot2))
    if size <= 50:
        return 4
    if size <= 100:
        return 10
    if size <= 200:
        return 12
    return 15

def distanceBucket(q): #1 angstrom bucketsize
    return (math.floor(q[0]),math.floor(q[1]),math.floor(q[2]),math.floor(q[3]))
def triangleBucket(q):
    return match.floor(q[4]/2)

def make_buckets(pdict):
    buckets = [[[[[]]]]]
    for p in pdict.getKeys():
        for q in pdict[p]:
            (i1,i2,i3,i4) = distanceBucket(q)
            i5 = triangleBucket(q)
            buckets[i1][i2][i3][i4][i5] = p
