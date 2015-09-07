import csv
from matplotlib import pyplot
f = open('../logs/EigenLog','r')
lol = list(csv.reader(f, delimiter='\t'))
i = 0
listByProteins = [[],[]]
Percentiles = [0,0,0,0,0,0,0,0,0,0]
for a in lol:
    if(len(a) <= 1 and len(listByProteins[i]) > 0):
        i = i + 1
        listByProteins.append([])
    elif(a[0] != ' EigenValues '):
        listByProteins[i] = listByProteins[i] + a
listByProteins.pop(0)
listByProteins.reverse()
listByProteins.pop(0)
listByProteins.reverse()
finalList = []
for a in listByProteins:
    a = [float(i) for i in a]
    m = max(a)
    tempList = [float(x)/m for x in a]
    for i in tempList:
        if(i <= .1):
            Percentiles[0] = Percentiles[0] + 1
        elif(i <= .2):
            Percentiles[1] = Percentiles[1] + 1
        elif(i <= .3):
            Percentiles[2] = Percentiles[2] + 1
        elif(i <= .4):
            Percentiles[3] = Percentiles[3] + 1
        elif(i<=.5):
            Percentiles[4] = Percentiles[4] + 1
	elif(i<=.6):
	    Percentiles[5] = Percentiles[5] + 1
	elif(i<=.7):
	    Percentiles[6] = Percentiles[6] + 1
	elif(i<=.8):
	    Percentiles[7] = Percentiles[7] + 1
	elif(i<=.9):
	    Percentiles[8] = Percentiles[8] + 1
	else:
	    Percentiles[9] = Percentiles[9] + 1
s = sum(Percentiles)
Percentiles = [float(i)/float(s) for i in Percentiles]
print('The current Percentages for each 10% grouping are: ')
print(Percentiles)
pyplot.bar([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],Percentiles,width = 0.8,bottom=0,hold = 'none')
