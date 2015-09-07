/*
dpAlign.c

Authors: Allen Holder, Jonathon Strouser, and Jonathan Taylor
Date: 6/2013

The DP portions of this code are adaptations of J. Strouser's code to 
align two sequences with affine gaping dynamic programming.  The original 
code was written as a mex file to work with Matlab, but this alteration
compiles with swig as an extension to python.

The assessment function is a re-write of J. Taylor's matlab code to calculate
a ROC curve and its area after a database alignment is undertaken.  This to
has been adapted to compile with swig as an extension to python.

Install:
	make; make install
   or
	swig -python dpAlign.i
	gcc -c -fPIC dpAlign.c dpAlign_wrap.c -I/usr/include/python2.7
	ld -shared dpAlign.o dpAlign_wrap.o -o _dpAlign.so
	cp ./dpAlign.py ..
        cp ./_dpAlign.so ..

Use in a python file (Example):
	import dpAlign
	# Create arrays for swig's wrapper
        seq1 = dpAling.doubleArray(lenSeq1)
	seq2 = dpAlign.doubleArray(lenSeq2)
        # populate the sequences with some numerical values, e.g. the
        # eigenvalues of the two encoded sequences.  Then the following
	# command returns a (swig)pointer to the alignment structure
        ptrAlignStruct = dpAlign.CalcGlobalAlignment(s1, lenS1, s2, lenS2, \
                                             gapOpenPen, gepExtendPen)
	# We can get the alignment infomation and free the structure's
	# memory with
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

        # Handling the assessment code is similar.  See EIGAs.py and
	# EIGAs_utils.py
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define NEGINF -1E6
#define MAX(x,y) (x) > (y) ? (x) : (y)
#include "dpAlign.h"

void affine_function(int input, double gapOpen, 
                     double gapExtend, double* temp){
   *temp = (double)((double)gapOpen + (double)(gapExtend*(double)input));
}

// Find location of the maximum of 3 values
int maxVal(double v1, double v2, double v3){
   double m;
   m = (MAX((MAX(v1, v2)), v3));
   if(m==v1)
      return 1;
   else if(m==v2)
      return 2;
   else if(m==v3)
      return 3;
   else{
      return -1;
   }
}

// Calculate an optimal, global alignment with affine gaping DP.
// Return and pointer to an alignInfo structure
alignInfo *CalcGlobalAlignment(double *seq1, int lenSeq1, double *seq2, 
               int lenSeq2, double gapOpen, double gapExtend) {
    
   double **matchScore; // pairwise values of matching residues
   int i, j; // indices
   int **index; // Keeps track of DPs trace through the matrix
   double **score, **gh, **gv; // DP scoring matrices that store the DP process
   double val1, val2;  // interim calculated values
   double temp; // temporary calculation
   double optValOfAlignment; // Returned value of optimal global alignment
   int *I1, *I2; // arrays to store the consecutive matched pairs
   int lenI1, lenI2; // size of the resulting matched pairs
   alignInfo *alignment=malloc(sizeof(alignInfo)); // structure that we return
   
   // Allocate space
   matchScore = (double **) calloc(lenSeq1+1, sizeof(double *));
   score = (double **) calloc(lenSeq1+1, sizeof(double *));
   gh = (double **) calloc(lenSeq1+1, sizeof(double *));
   gv = (double **) calloc(lenSeq1+1, sizeof(double *));
   index = (int **) calloc(lenSeq1+1, sizeof(int *));
   I1 =  (int *) calloc(1000,sizeof(int));
   I2 =  (int *) calloc(1000,sizeof(int));
    
   for(i=0; i<lenSeq1+1; i++){
      matchScore[i] = (double *) calloc(lenSeq2+1, sizeof(double));
      score[i] = (double *) calloc(lenSeq2+1, sizeof(double));
      gh[i] = (double *) calloc(lenSeq2+1, sizeof(double));
      gv[i] = (double *) calloc(lenSeq2+1, sizeof(double));
      index[i] = (int *) calloc(lenSeq2+1, sizeof(int));
      if(index[i] == NULL){
          printf("OUT OF MEMORY\n");
          exit(1);
       }
   }
    
   // Calcluate the EIGAs score matrix for the two sequences
   for(i=0; i<lenSeq1; i++)
      for(j=0; j<lenSeq2; j++) {
         matchScore[i][j] = -fabs(seq1[i]-seq2[j]) / (seq1[i] + seq2[j]);
      } 
    
   // Initialize the boarders of the DP process
   score[0][0] = gh[0][0] = gv[0][0] = 0.0;
   for(j=1; j < lenSeq2+1; j++){
      score[0][j] = 0.0-9999.0;
      affine_function(j, gapOpen, gapExtend, &temp);
      gv[0][j] = 0.0 - temp;
      gh[0][j] = 0.0-9999.0;
      index[0][j] = 1;
   }
   for(i=1; i < lenSeq1+1; i++){
      score[i][0] = 0.0-9999.0;
      affine_function(i, gapOpen, gapExtend, &temp);
      gh[i][0] = 0.0 - temp;
      gv[i][0] = 0.0-9999.0;
      index[i][0] = 2;
   }
    
   // Forward sweep of the DP recursion
   for(i=1; i < lenSeq1+1; i++)
      for(j=1; j < lenSeq2+1; j++){
         val1 = gv[i][j-1] - gapExtend;
         val2 = gh[i-1][j] - gapExtend;
         score[i][j] = matchScore[i-1][j-1] 
                       + (MAX((MAX(score[i-1][j-1], gv[i-1][j-1])),
                               gh[i-1][j-1]));
         gv[i][j] = (MAX((MAX((score[i][j-1] - (gapOpen+gapExtend)), val1)), 
                     gh[i][j-1] - (gapOpen+gapExtend)));
         gh[i][j] = (MAX((MAX((score[i-1][j] - (gapOpen+gapExtend)), 
                     (gv[i-1][j] - (gapOpen+gapExtend)))), val2));
         index[i][j] = maxVal(gv[i][j], gh[i][j], score[i][j]);
      }
    
   // Store the global alignment value before freeing memory
   optValOfAlignment = (MAX((MAX(score[lenSeq1][lenSeq2],
                        gh[lenSeq1][lenSeq2])), gv[lenSeq1][lenSeq2]));

   // Calculate the matched pairs
   i = lenSeq1;
   j = lenSeq2;
   lenI1 = lenI2 = 0;
   while((i>0) || (j>0)){
      if(index[i][j] == 3){
         I1[lenI1++] = i;
         I2[lenI2++] = j;
         i--;j--;
      }
      else if(index[i][j] == 1){
         j--;
      }
      else {
        i--;
      }
   }
   alignment->DPscore = optValOfAlignment;
   alignment->lenSeq1 = lenI1;
   alignment->lenSeq2 = lenI2;
   alignment->alignSeq1 = (int *) malloc(lenI1*sizeof(int));
   alignment->alignSeq2 = (int *) malloc(lenI2*sizeof(int));
   for(i=0; i<lenI1; i++) {
      alignment->alignSeq1[i] = I1[lenI1-i-1];
   }
   for(i=0; i<lenI2; i++) {
      alignment->alignSeq2[i] = I2[lenI2-i-1];
   }
   
   // Free the memory
   for(i = 0; i < lenSeq1+1; i++){
     free(matchScore[i]);
     free(score[i]);
     free(gh[i]);
     free(gv[i]);
     free(index[i]);
   }
   free(matchScore);
   free(score);
   free(gh);
   free(gv);
   free(index);
   free(I1);
   free(I2);

   // We are done
   return(alignment);
}

void freeAlignInfo(alignInfo *alignment) {
   free(alignment->alignSeq1);
   free(alignment->alignSeq2);
   free(alignment);
}

double getDPVal(alignInfo *alignment) {
   return(alignment->DPscore);
}

int getSeqLen(alignInfo *alignment, int seqNum) {
   if (seqNum == 1) {
      return(alignment->lenSeq1);
   } else if (seqNum == 2) {
      return(alignment->lenSeq2);
   } else {
      return(0);
   }
}

int getSeq(alignInfo *alignment, int seqNum, int subscript) {
   if (seqNum == 1) {
      if ((subscript >= 0) && (subscript < alignment->lenSeq1)) {
         return(alignment->alignSeq1[subscript]);
      } else {
         return(0);
      }
   } else if (seqNum == 2) {
      if ((subscript >= 0) && (subscript < alignment->lenSeq2)) {
         return(alignment->alignSeq2[subscript]);
      } else {
         return(0);
      }
   } else {
      return(0);
   }
}

// Assess an alignment to find family agreement and it's ROC curve.
// Returns a pointer to structure assessInfo
assessInfo *assessAlignment(int *familyID, double *alignScore, 
                               int numProt) {

   int numCorrect; // number of families identified by best alignment score
   double maxScore, minScore; // max and min alignment scores
   double bstScore; // best score for each protein
   int bstProt; // closest associated protein for each protein
   int prt1, prt2, b; // iterators
   int delta; // refinement over which we calculate the ROC curve
   int sizeRefinement; // number of points in the partition of the ROC curve
   int **results; // stores counts for ROC the curve
   double *TPR, *FPR; // True and false positive rates for the ROC curve
   double area; // area under the ROC curve
   double dist, bestDist, bestDiscrn, bestTPR, bestFPR; // used to calculate
                                                        // best point on ROC
   assessInfo *assessment=malloc(sizeof(assessInfo)); // structure we return

   // Calculate how many nearest pairs agree on family ID; along the way
   //find the largest and smallest values
   numCorrect = 0;
   maxScore = -100000000;
   minScore = 100000000;
   for (prt1 = 0; prt1 < numProt; prt1++) {
      bstScore = -100000000;
      for (prt2 = 0; prt2 < numProt; prt2++) {
         if (prt1 != prt2) {
            if (alignScore[prt1*numProt+prt2] > bstScore) {
               bstScore = alignScore[prt1*numProt+prt2];
               bstProt = familyID[prt2];
            }
            if (alignScore[prt1*numProt+prt2] > maxScore) {
               maxScore = alignScore[prt1*numProt+prt2];
            }
            if (alignScore[prt1*numProt+prt2] < minScore) {
               minScore = alignScore[prt1*numProt+prt2];
            }
         }
      }
      if (familyID[prt1] == bstProt) {
         numCorrect = numCorrect + 1;
      }
   }

   // Calculate a ROC curve and its area
   //
   // Allocate space
   delta = 1000;
   sizeRefinement = floor((maxScore-minScore)*delta+1);
   results = (int **) calloc(sizeRefinement, sizeof(int *));
   TPR = (double *) calloc(sizeRefinement, sizeof(double));
   FPR = (double *) calloc(sizeRefinement, sizeof(double));
   for (b = 0; b < sizeRefinement; b++) {
      results[b] = (int *) calloc(4, sizeof(int));
      if (results[b] == NULL) {
         printf("OUT OF MEMORY\n");
         exit(1);
      }
   }
   // Initialize results to zero
   for (b = 0; b < sizeRefinement; b++) {
      results[b][0] = 0;
      results[b][1] = 0;
      results[b][2] = 0;
      results[b][3] = 0;
   }
   // Calculate the true/false positives and the true/false negatives
   for (b = 0; b < sizeRefinement; b++) {
      for (prt1 = 0; prt1 < numProt; prt1++) {
         for (prt2 = 0; prt2 < prt1; prt2++) {
            if (alignScore[prt1*numProt+prt2] >= minScore+(1.0*b/delta)) {
               if (familyID[prt1] == familyID[prt2]) {
                  results[b][0] = results[b][0]+1;
               } else if (familyID[prt1] != familyID[prt2]) {
                  results[b][1] = results[b][1]+1;
               } else {
                  printf("ERROR calculating ROC curve (1)\n");
               } 
            } else if (alignScore[prt1*numProt+prt2] < minScore+(1.0*b/delta)) {
               if (familyID[prt1] == familyID[prt2]) {
                  results[b][2] = results[b][2]+1;
               } else if (familyID[prt1] != familyID[prt2]) {
                  results[b][3] = results[b][3]+1;
               } else {
                  printf("ERROR calculating ROC curve (2)\n");
               }
            } else {
               printf("ERROR calculating ROC curve (3)\n");
            }
         }
      }
   } 

   // Calculate the true and false positive rates:
   bestDist = 10;
   bestDiscrn = 0;
   for (b = 0; b < sizeRefinement; b++) {
      TPR[b] = 1.0*results[b][0]/(results[b][0] + results[b][2]);
      FPR[b] = 1.0*results[b][1]/(results[b][1] + results[b][3]);
      dist = sqrt((1-TPR[b])*(1-TPR[b]) + FPR[b]*FPR[b]);
      if (dist < bestDist) {
         bestDist = dist;
         bestDiscrn = minScore+(1.0*b/delta);
         bestTPR = TPR[b];
         bestFPR = FPR[b];
      }
   }
   // Integrate
   area = 0;
   for (b = 1; b < sizeRefinement; b++) {
      area = area + TPR[b]*(FPR[b-1]-FPR[b]);
   }

   // populate the return structure
   assessment->numCorrect = numCorrect;
   assessment->partitionSize = sizeRefinement;
   assessment->ROCarea = area;
   assessment->maxScore = maxScore;
   assessment->minScore = minScore;
   assessment->bestDiscernmentVal = bestDiscrn;
   assessment->bestTPR = bestTPR;
   assessment->bestFPR = bestFPR;
   assessment->TPR = (double *) malloc(sizeRefinement*sizeof(double));
   assessment->FPR = (double *) malloc(sizeRefinement*sizeof(double));
   for (b = 0; b < sizeRefinement; b++) {
      assessment->TPR[b] = TPR[b];
      assessment->FPR[b] = FPR[b];
   }

   // Free the memory
   for (b = 0; b < sizeRefinement; b++) {
      free(results[b]);
   }
   free(results);
   free(TPR);
   free(FPR);

   // return the assessment structure
   return(assessment);
}

double getTPR(assessInfo *assessment, int i) {
   return(assessment->TPR[i]);
}

double getFPR(assessInfo *assessment, int i) {
   return(assessment->FPR[i]);
}

int getNumCorrect(assessInfo *assessment) {
   return(assessment->numCorrect);
}

int getPartitionSize(assessInfo *assessment) {
   return(assessment->partitionSize);
}

double getROCarea(assessInfo *assessment) {
   return(assessment->ROCarea);
}

double getMaxScore(assessInfo *assessment) {
   return(assessment->maxScore);
}

double getMinScore(assessInfo *assessment) {
   return(assessment->minScore);
}

double getBestDiscernmentVal(assessInfo *assessment) {
   return(assessment->bestDiscernmentVal);
}

double getBestTPR(assessInfo *assessment) {
   return(assessment->bestTPR);
}

double getBestFPR(assessInfo *assessment) {
   return(assessment->bestFPR);
}

void freeAssessInfo(assessInfo *assessment){
   free(assessment->TPR);
   free(assessment->FPR);
   free(assessment);
}
