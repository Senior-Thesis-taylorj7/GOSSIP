#!/usr/bin/python

#
# EIGAs, ver 0.6
#
# Author: Allen Holder
# Date: 6-7/2013
#

# This is the primary gateway to EIGAs.
# Set the following appropriately and then call this python script.

from EIGAs_main import *

# Settings

#dataBaseDir = '../data/Proteus300/' # Directory that holds the pdb files
#dataBaseDir = '../data/test/' # Directory that holds the pdb files
#parsedDataFile = '../tmp/test.eigas' # A pyhon cPickle file name:
dataBaseDir = '../data/Skolnick40/' # Directory that holds the pdb files
parsedDataFile = '../tmp/Skolnick40.eigas' # A pyhon cPickle file name:
#parsedDataFile = '../tmp/Proteus300.eigas' # A pyhon cPickle file name:
                                           # can be read from if data has
                                           # already been parsed, or,
                                           # will be written two if encoding
                                           # from the pdb files
parseData = True # Set to True if you want to read the coordinates from
                 # the pdf/ent files and encode with these newly read 
                 # coordinates
saveParsedData = False # Set to True if you want to save the parsed data, i.e.
                       # a list of all lambda sequences for their respective
                       # chains
kappa = 17 # kappa value, cuttoff for smooth contact to be zero
gapOpen = 0.9 # gap open penalty for DP
gapExtend = 0.4 # gap extension penalty for DP
numProcessors = 20 # Number of processors to distribute the pairswise
                  # comparisons over
classificationFile = '../data/SCOP/dir.cla.scop.txt_1.75' # The SCOP
                                                          # classification file

# Send to EIGAs.  Returns an assessment dictionary
print('\n\n***** EIGAs *****\n')
assessment = EIGAs(dataBaseDir, parsedDataFile, parseData, saveParsedData, \
         kappa, gapOpen, gapExtend, numProcessors, classificationFile)

# Print some results
print('EIGAs -> '+str(assessment['numCorrect'])+'/'+\
      str(assessment['numProteins'])+\
      ' were correct as'+' measured by best score')
print('EIGAs -> ROC curve area was '+str(assessment['ROCarea']))
print('EIGAs -> minScore <= bestDiscernment <= maxScore were ')
print('\t'+str(assessment['minScore'])+' <= '+\
      str(assessment['bestDiscernmentVal'])+' <= '+\
      str(assessment['maxScore']))
print('EIGAs -> Best TPR and FPR were '+str(assessment['bestTPR'])+\
         ' and '+str(assessment['bestFPR']))
dist = sqrt((1-assessment['bestTPR'])**2+assessment['bestFPR']**2)
print('\t distance to perfect (0,1) on ROC curve was '+str(dist))
print('\t ROC graph is stored as ../tmp/ROC.png')

print('\n\n***** EIGAs Done *****')
