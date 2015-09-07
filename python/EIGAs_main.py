# Import python extensions for EIGAs
import os
import EIGAs_utils
import cPickle
import itertools
import time
#import matplotlib
from pylab import *
from multiprocessing import Pool
   
def EIGAs(dataBaseDir, parsedDataFile, parseData, saveParsedData, \
               kappa, gapOpen, gapExtend, numProcessors, classificationFile):

   # Open Log File
   if not os.path.exists('../logs/'):
      os.makedirs('../logs/')
   logFile = open('../logs/runlog','w')
   HistFile = open('../logs/EigenLog','w')
   print('EIGAs -> Run details in '+'../logs/runlog')
   
   # Parse Proteins
   print('EIGAs -> Collecting and Encoding the Proteins')
   ticParse = time.time()
   if parseData:
      logFile.write('Processing files from '+dataBaseDir+'\n')
      tic = time.time()
      Prot = EIGAs_utils.readProteins(logFile, dataBaseDir, classificationFile)
      logFile.write('\t required '+str(time.time()-tic)+' seconds\n')
   
      # Encode each chain as a list of eigenvalues
      tic = time.time()
      Prot = EIGAs_utils.encode(Prot, kappa, logFile, HistFile)
      logFile.write('\t required '+str(time.time()-tic)+' seconds\n')
   
      # Save the datafile
      if saveParsedData:
         logFile.write('Saving parsed data in '+parsedDataFile+'\n')
         tic = time.time()
         pklFile = open(parsedDataFile, 'wb')
         cPickle.dump(Prot, pklFile)
         pklFile.close()
         logFile.write('\t required '+str(time.time()-tic)+' seconds\n')
   else:
      # Load the datafile
      if os.path.exists(parsedDataFile):
         logFile.write('Loading parsed data from '+parsedDataFile+'\n')
         tic = time.time()
         pklFile = open(parsedDataFile, 'rb')
         Prot = cPickle.load(pklFile)
         pklFile.close()
         logFile.write('\t required '+str(time.time()-tic)+' seconds\n')
      else:
         logFile.write('No data file found\n')
   print('\t\t required '+str(time.time()-ticParse)+' seconds')
   
   # Perform pairwise alignments
   print('EIGAs -> Conducting pairwise alignments')
   logFile.write('\n ***** EIGAs Beginning Pairwise Alignments *****\n')
   ticAlign = time.time()
   if numProcessors > 1:
      pool = Pool(processes=numProcessors)
   alignment = {}
   for (prtName1, prtName2) in itertools.combinations(Prot, 2):
      # Check to make sure we have legitimate chains
      if Prot[prtName1].model[0].chain[0].legitChain and \
         Prot[prtName2].model[0].chain[0].legitChain:
         # Send off to align with DP
         if numProcessors > 1:
            alignment[(prtName1, prtName2, 'score')], \
               alignment[(prtName1, prtName2, 'seq1')], \
               alignment[(prtName1, prtName2, 'lenSeq1')], \
               alignment[(prtName1, prtName2, 'seq2')], \
               alignment[(prtName1, prtName2, 'lenSeq2')] = \
               pool.apply_async(EIGAs_utils.align, \
               (Prot[prtName1].model[0].chain[0].lambdaVal, \
                Prot[prtName2].model[0].chain[0].lambdaVal, \
                gapOpen, gapExtend) ).get()
         elif  numProcessors == 1:
            alignment[(prtName1, prtName2, 'score')], \
               alignment[(prtName1, prtName2, 'seq1')], \
               alignment[(prtName1, prtName2, 'lenSeq1')], \
               alignment[(prtName1, prtName2, 'seq2')], \
               alignment[(prtName1, prtName2, 'lenSeq2')] = \
               EIGAs_utils.align(\
                  Prot[prtName1].model[0].chain[0].lambdaVal, \
                  Prot[prtName2].model[0].chain[0].lambdaVal, \
                  gapOpen, gapExtend)
         else:
            logFile.write('Error, numProcessors is unknown\n')
            return {}
         # normalize the score
         alignment[(prtName1, prtName2, 'score')] = \
            alignment[(prtName1, prtName2, 'score')] / (0.5 * \
            len(Prot[prtName1].model[0].chain[0].lambdaVal) + \
            len(Prot[prtName2].model[0].chain[0].lambdaVal) )
         # Calculate the sCMO score
         alignment[(prtName1, prtName2, 'sCMO')] = \
            EIGAs_utils.calc_sCMO(\
               Prot[prtName1].model[0].chain[0].contactMat, \
               alignment[(prtName1, prtName2, 'seq1')],\
               Prot[prtName2].model[0].chain[0].contactMat,\
               alignment[(prtName1, prtName2, 'seq2')])
         # store the symmetric comparison
         alignment[(prtName2, prtName1, 'score')] = \
            alignment[(prtName1, prtName2, 'score')]
         alignment[(prtName2, prtName1, 'seq1')] = \
            alignment[(prtName1, prtName2, 'seq2')]
         alignment[(prtName2, prtName1, 'lenSeq1')] = \
            alignment[(prtName1, prtName2, 'lenSeq2')]
         alignment[(prtName2, prtName1, 'seq2')] = \
            alignment[(prtName1, prtName2, 'seq1')]
         alignment[(prtName2, prtName1, 'lenSeq2')] = \
            alignment[(prtName1, prtName2, 'lenSeq1')]
         alignment[(prtName2, prtName1, 'sCMO')] = \
            alignment[(prtName1, prtName2, 'sCMO')]
         # log the pairwise comparison
         logFile.write(str(prtName1)+', '+str(prtName2)+' has value '+\
                       str(alignment[(prtName1, prtName2, 'score')])+'\n')
   
   # Close the multiprocessing pool
   if numProcessors > 1:
      pool.close()
      pool.join()
   logFile.write('\t required '+str(time.time()-ticAlign)+' seconds over '\
         +str(numProcessors)+' processors\n')
   logFile.write('\n ***** EIGAs Done With Pairwise Alignments *****\n')
   print('\t\t required '+str(time.time()-ticAlign)+' seconds')
   
   # Assess the alignments
   print('EIGAs -> Assessing the Alignments')
   logFile.write('\nAssessing the alignments\n')
   ticAssess = time.time()
   assessment={}
   assessment['numCorrect'], assessment['TPR'], assessment['FPR'], \
      assessment['ROCarea'], assessment['partitionSize'], \
      assessment['maxScore'], assessment['minScore'], \
      assessment['bestDiscernmentVal'], assessment['bestTPR'], \
      assessment['bestFPR'] = EIGAs_utils.assessAlignment(Prot, alignment)
   assessment['numProteins'] = len(Prot)
   
   logFile.write('\t required '+str(time.time()-tic)+' seconds\n')
   logFile.write('\t '+str(assessment['numCorrect'])+'/'+str(len(Prot))+\
      ' were correct as measured by best score\n')
   logFile.write('\t ROC curve area was '+str(assessment['ROCarea'])+'\n')
   logFile.write('\t minScore <= bestDiscernment <= maxScore were \n')
   logFile.write('\t\t'+str(assessment['minScore'])+' <= '+\
         str(assessment['bestDiscernmentVal'])+' <= '+\
         str(assessment['maxScore'])+'\n')
   logFile.write('\t best TPR and FPR were '+str(assessment['bestTPR'])+\
         ' and '+str(assessment['bestFPR'])+'\n')
   dist = sqrt((1-assessment['bestTPR'])**2+assessment['bestFPR']**2)
   logFile.write('\t\t distance to perfect (0,1) on ROC curve was '+\
                 str(dist)+'\n')
   print('\t\t required '+str(time.time()-ticAssess)+' seconds')
   
   # Close Log File
   logFile.close()
   HistFile.close()

   return assessment
