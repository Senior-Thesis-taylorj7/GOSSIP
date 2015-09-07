typedef struct { 
   double DPscore;
   int lenSeq1, lenSeq2;
   int *alignSeq1, *alignSeq2;
} alignInfo;

typedef struct {
   int numCorrect;
   double *TPR, *FPR;
   double ROCarea;
   int partitionSize;
   double maxScore, minScore;
   double bestDiscernmentVal, bestTPR, bestFPR;
} assessInfo;

void affine_function(int input, double gapOpen,
                     double gapExtend, double* temp);

int maxVal(double v1, double v2, double v3);

alignInfo *CalcGlobalAlignment(double *seq1, int lenSeq1, double *seq2,
               int lenSeq2, double gapOpen, double gapExtend);

void freeAlignInfo(alignInfo *alignment);

double getDPVal(alignInfo *alignment);

int getSeqLen(alignInfo *alignment, int seqNum);

int getSeq(alignInfo *alignment, int seqNum, int subscript);

assessInfo *assessAlignment(int *familyID, double *alignScore, int numProt);

double getTPR(assessInfo *assessment, int i);

double getFPR(assessInfo *assessment, int i);

int getNumCorrect(assessInfo *assessment);

int getPartitionSize(assessInfo *assessment);

double getROCarea(assessInfo *assessment);

double getMaxScore(assessInfo *assessment);

double getMinScore(assessInfo *assessment);

double getBestDiscernmentVal(assessInfo *assessment);

double getBestTPR(assessInfo *assessment);

double getBestFPR(assessInfo *assessment);

void freeAssessInfo(assessInfo *assessment);
