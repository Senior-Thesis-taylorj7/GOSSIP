#
# EIGAs_utils
#
# Author: Allen Holder
# Date: 6/2013
#

import os
import itertools
import Bio.PDB
import warnings
from numpy import *
import matplotlib
from pylab import *
import dpAlign

#
# readProtiens
# 
# reads protein files in the data directory and creates EIGAs structure, Prot
#
def readProteins(logFile, dataBaseDir, classificationFile):
   if os.path.exists(dataBaseDir):
      fileList = os.listdir(dataBaseDir)
   else:
      logFile.write("Input data path does not exist.\n")
   Prot = {}
   # Read the classification file
   classFile = open(classificationFile, 'r')
   classification = classFile.read()
   classFile.close()
   # Temporarily turn off warnings since biopython gives many due to
   # discontinuous chains
   with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      # Set the parser with permission to ignore file errors
      parser = Bio.PDB.PDBParser(PERMISSIVE=1)
      # Begin parsing
      for fname in fileList:
         # Only grab the PDB files ending with pdf or ent.
         # Create a list of filenames without extensions.
         if fname.endswith('.pdb') or fname.endswith('.ent'):
            prtName = fname[0:len(fname)-4]
            # Use biopython to parse the pdb file
            logFile.write('Processing '+fname+'\n')
            # Create a protein structure
            Prot[prtName] = parser.get_structure(prtName,dataBaseDir+fname)
            # Get the protein family
            tmpLine = classification[classification.find(str(prtName).lower()):]
            Prot[prtName].familyID = \
                  tmpLine[tmpLine.find('fa=')+3:tmpLine.find('dm=')-1]
            # store the models for the structure
            Prot[prtName].model = [model for model in Prot[prtName]]
            for modNum in range(len(Prot[prtName].model)):
               logFile.write('\tModel '+str(modNum)+'\n')
               # store the chains for each model
               Prot[prtName].model[modNum].chain = \
                     [chain for chain in Prot[prtName].model[modNum]]
               # collect the pertinant information for each chain
               for chainNum in range(len(Prot[prtName].model[modNum].chain)):
                  logFile.write('\t\tChain '+str(chainNum)+', ')
                  # Count the C-alpha atoms
                  numRes = 0
                  for residue in Prot[prtName].model[modNum].chain[chainNum]:
                     if residue.has_id('CA'):
                        numRes = numRes + 1
                  logFile.write('\t\t\t# Residues = '+str(numRes)+'\n')
                  # Read the coordinates and bfactors in local numpy arrays
                  coords = zeros((numRes,3))
                  bfactor = zeros((numRes,1))
                  i = 0
                  for residue in Prot[prtName].model[modNum].chain[chainNum]:
                     if residue.has_id('CA'):
                        (coords[i,0], coords[i,1], coords[i,2]) = \
                              residue['CA'].get_coord()
                        bfactor[i] = residue['CA'].get_bfactor()
                        i = i + 1
                  # Store them in the protein container
                  Prot[prtName].model[modNum].chain[chainNum].coords = coords
                  Prot[prtName].model[modNum].chain[chainNum].bfactor = bfactor
         else:
            logFile.write('Protein '+fname+' is not pdb or ent\n')

   logFile.write('***** EIGAs Done Parsing *****\n')

   # Return the collected information
   return Prot

#
# encode
#
# Calculate the smooth contact matrices for each chain with respect to a
# cutoff value kappa, and then encode the resdiues as eigenvalues
#
def encode(Prot, kappa, logFile, HistFile):
   logFile.write('\n\n ***** EIGAs Encoding *****\n')
   HistFile.write('Data')
   for prtName in Prot:
      logFile.write('Encoding '+prtName+' ')
      HistFile.write('Protein : '+prtName+'\n EigenValues \n')
      for modNum in range(len(Prot[prtName].model)):
         for chainNum in range(len(Prot[prtName].model[modNum].chain)):
            coords = Prot[prtName].model[modNum].chain[chainNum].coords
            numRes,numCoords = coords.shape
            legitChain = True
            if numCoords == 3 and numRes > 0:
               C = zeros((numRes, numRes))
               for i in xrange(numRes):
                  for j in xrange(i+1):
                     dist = linalg.norm(coords[i,:] - coords[j,:])
                     if dist <= kappa and dist >= 0:
                        C[i,j] = 1 - dist/kappa
                     elif dist > kappa:
                        C[i,j] = 0
                     else:
                        logFile.write('\tUnknown distance between residues\n')
            elif numCoords != 3:
               logFile.write('\tWrong number of coordinates, skipping\n')
               legitChain = False
            else:
               logFile.write('\tNo residues, skipping\n')
               legitChain = False
            # encode the residues if we have a legitimate chain
            if legitChain:
               # Factor the smooth contact map
               D,U = linalg.eigh(C) #This calculates the Eigenvectors and Eigenvalues of the relevant contact map
               # Project the residue onto the nearest eigenspace in terms of
               # contact angle
               lambdaIndx = abs(dot(U, \
                     dot(sign(diag(D)), sqrt(abs(diag(D)))))).argmax(axis=1)
               Prot[prtName].model[modNum].chain[chainNum].lambdaVal \
                     = D[lambdaIndx]
               D[lambdaIndx].tofile(HistFile,'\t','%.10s')
               HistFile.write('\n')
	       Prot[prtName].model[modNum].chain[chainNum].contactMat = C
            Prot[prtName].model[modNum].chain[chainNum].legitChain \
                  = legitChain
      logFile.write('\n')

   # Done encoding
   logFile.write('\n ***** EIGAs Done Encoding *****\n')
   return Prot

#
# align
#
# Calculate the pairwise score matrix of two lambda sequences and align with
# biopythons dynamic programming
#
def align(s1, s2, gapOpen, gapExtend):
   # instantiate the data type to hand to the wraped C code
   t1 = dpAlign.doubleArray(len(s1))
   t2 = dpAlign.doubleArray(len(s2))
   # populate the sequences with the correct values
   for i in range(len(s1)):
      t1[i] = s1[i]
   for i in range(len(s2)):
      t2[i] = s2[i]
   # Calculate the alignment with DP
   alignStruct = dpAlign.CalcGlobalAlignment(t1, len(s1), t2, len(s2), \
                                         gapOpen, gapExtend)
   score = dpAlign.getDPVal(alignStruct)
   lenSeq1 = dpAlign.getSeqLen(alignStruct, 1)
   lenSeq2 = dpAlign.getSeqLen(alignStruct, 2)
   seq1 = {}
   seq2 = {}
   for i in range(lenSeq1):
      seq1[i] = dpAlign.getSeq(alignStruct, 1, i)
   for i in range(lenSeq2):
      seq2[i] = dpAlign.getSeq(alignStruct, 2, i)
   dpAlign.freeAlignInfo(alignStruct)
   return score, seq1, lenSeq1, seq2, lenSeq2

#
# calc_sCMO
#
def calc_sCMO(C1, s1, C2, s2):
   ind1 = zeros((len(s1)),int)
   ind2 = zeros((len(s2)),int)
   for i in range(len(s1)):
      ind1[i] = s1[i]-1
   for i in range(len(s2)):
      ind2[i] = s2[i]-1
   C1upper = triu(C1[ix_(ind1,ind1)],2)
   C2upper = triu(C2[ix_(ind2,ind2)],2)
   C1C2 = C1upper*C2upper
   C1plusC2 = C1upper.sum() + C2upper.sum()
   if C1plusC2 > 0:
      return C1C2/C1plusC2
   elif C1plusC2 == 0:
      return 0 # denotes that there is no contact
   else:
      return -2 # denotes an error

def assessAlignment(Prot, alignment):
   familyID = dpAlign.intArray(len(Prot))
   alignScore = dpAlign.doubleArray(len(Prot)**2)
   # build a dictionary of protein names to ensure consistency
   i = 0
   nameList = {}
   for prtName in Prot:
      nameList[i] = prtName
      i = i + 1
   # build the data structures for the wrapped C code
   for i in range(len(Prot)):
      familyID[i] = int(Prot[nameList[i]].familyID)
      for j in range(len(Prot)):
         if i != j:
            alignScore[i*len(Prot)+j] = \
               alignment[(nameList[i], nameList[j], 'score')]
         else:
            alignScore[i*len(Prot)+j] = 0
   assessmentStruct = dpAlign.assessAlignment(familyID, alignScore, len(Prot))
   numCorrect = dpAlign.getNumCorrect(assessmentStruct)
   ROCarea = dpAlign.getROCarea(assessmentStruct)
   partitionSize = dpAlign.getPartitionSize(assessmentStruct)
   TPR = zeros((partitionSize, 1))
   FPR = zeros((partitionSize, 1))
   for i in range(partitionSize):
      TPR[i] = dpAlign.getTPR(assessmentStruct, i)
      FPR[i] = dpAlign.getFPR(assessmentStruct, i)
   maxScore = dpAlign.getMaxScore(assessmentStruct)
   minScore = dpAlign.getMinScore(assessmentStruct)
   bestDiscernmentVal = dpAlign.getBestDiscernmentVal(assessmentStruct)
   bestTPR = dpAlign.getBestTPR(assessmentStruct)
   bestFPR = dpAlign.getBestFPR(assessmentStruct)
   dpAlign.freeAssessInfo(assessmentStruct)
   plot(FPR, TPR, linewidth=2.0)
   plot(arange(0.0,1.01,0.01), arange(0.0,1.01,0.01),\
        linewidth=1.0, linestyle=':')
   title('ROC Curve')
   axis([-0.1, 1.1, -0.1, 1.1])
   xlabel('False Positive Rate (FPR)')
   ylabel('True Positive Rate (TPR)')
   if not os.path.exists('../tmp/'):
      os.makedirs('../tmp/')
   savefig("../tmp/ROC.png")
   return numCorrect, TPR, FPR, ROCarea, partitionSize,\
          maxScore, minScore, bestDiscernmentVal, bestTPR, bestFPR
