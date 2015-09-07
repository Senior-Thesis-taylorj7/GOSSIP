%module dpAlign
%{
#include "dpAlign.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

extern alignInfo *CalcGlobalAlignment(double *seq1, int lenSeq1, double *seq2,  
               int lenSeq2, double gapOpen, double gapExtend);

extern void freeAlignInfo(alignInfo *alignment);

extern double getDPVal(alignInfo *alignment);

extern int getSeqLen(alignInfo *alignment, int seqNum);

extern int getSeq(alignInfo *alignment, int seqNum, int subscript);

extern assessInfo *assessAlignment(int *familyID, double *alignScore, 
                                   int numProt);

extern double getTPR(assessInfo *assessment, int i);

extern double getFPR(assessInfo *assessment, int i);

extern int getNumCorrect(assessInfo *assessment);

extern int getPartitionSize(assessInfo *assessment);

extern double getROCarea(assessInfo *assessment);

extern double getMaxScore(assessInfo *assessment);

extern double getMinScore(assessInfo *assessment);

extern double getBestDiscernmentVal(assessInfo *assessment);

extern double getBestTPR(assessInfo *assessment);

extern double getBestFPR(assessInfo *assessment);

extern void freeAssessInfo(assessInfo *assessment);
