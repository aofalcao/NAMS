#ifndef AWTS_H
#define AWTS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "nams.h"




int GenerateAtomWeights(char *fname1, char *fname2, Params *parms, bool J, int *filt_db, int nruns, char *wtsFname, char *iwtsFname, float spread); 
int WeightsEvaluator(float ***matrices, int **mat, float *wts, int nmols, int natoms1, int *natoms2, float *scores, bool J);



#ifdef __cplusplus
}
#endif



#endif