/*
This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published on the official site of Open Source Initiative
and attached above.

Copyright (C) 2013, Andre Falcao and Ana Teixeira, University of Lisbon - LaSIGE

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Please cite the authors in any work or product based on this material:

AL Teixeira, AO Falcao. 2013. A non-contiguous atom matching structural similarity function. J. Chem. Inf. Model. DOI: 10.1021/ci400324u.
*/

#ifndef NAMS_H
#define NAMS_H

#ifdef __cplusplus
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "hungarian.h"


#define MAX_LIN_SIZ 1024
#define ABA_ELEMS 10
#define MAX_LEVELS 128
#define NAMS_VERSION "0.99.20170104B"

typedef struct aba {
	int numb1, nrings1, chir1, numb2, nrings2, chir2, inring, arom, order, dbcistrans;
} AbaBond;

typedef struct theparams{
    float BS_ALPHA;
    float ANRINGS_FAC;	   //  #number of rings an atom belongs to
    float ACHIR_FAC;		//  #chiral atom
    float DBSTEREO_FAC;	    //  #double bond stereo
    float BRING_FAC;		//  #bond in ring
    float BAROM_FAC;		//  #bond aromaticity
    float BORDER_FAC;		//  #bond order
    float PEN;
    int *LU_ELEMS;		// a lookup table for discovering the indices of the elements from the matrix from the atomic numbers
    float **ELEMS_DISTS;
}Params;

typedef struct mol {
	char mname[32];
	int  cid;
	int  molwt;					//molecular weight
	int  natoms;
	int  nbonds;
	int  naba_types;
	int      *nlevels;          //Number of levels in each bond (one value per atom)
	AbaBond *aba_types;		//the ababonds in themselves
	float    *weights;          //the weight to aplly to each atom - should most of the times be a vector of ones	
	int *abas;		            //these are the indices of the AbaBonds for each atom. They must match exactly with the levels
	int *levels;                //vector with levels for each ababond for each atom
} MolInfo;


//this is a plug in struct for sorting similarity between bonds of 2 molecules
typedef struct simbnd {
	int sim, bnd1, bnd2;
}SimBond;

bool cid_in_flist(int cid, int *filt_list);
int setElemsDists(Params *parms, int adm);
int setParams(Params *parms, int ADM, float BS_ALPHA,  float ANRINGS_FAC,  float ACHIR_FAC, float DBSTEREO_FAC, float BRING_FAC, float BAROM_FAC, float BORDER_FAC, float PEN);
int readMolInfo( FILE *fil, MolInfo *mi); 
int readMolInfoBin( FILE *fil, MolInfo *mi); 
int blankReadMolFile( FILE *fil); 
int blankReadMolFileBin( FILE *fil); 
int showMolInfo(MolInfo *mi);
int readSimAtomMAtrix(char *fname); 
int compare_aba_bonds(AbaBond *aba1, AbaBond *aba2, Params *parms);
int *getAbaBondsCompMatrix(MolInfo *mi1, MolInfo *mi2, Params *parms);
int *calcBondLevelsMatrix(Params *parms);
int simbond_compare( const void *sim1, const void *sim2 );
int matchBonds(int a1, int a2, MolInfo *mi1, MolInfo *mi2, int *aba_ms, int *blev_mat); 
int *calcAtomMatchingMatrix(MolInfo *mi1, MolInfo *mi2, int *aba_match_scores, int *blev_mat, float pen);
int nams_runner(MolInfo *mi1, MolInfo *mi2, Params *parms, int *blev_mat, bool M, bool A, float *wts);
int freeMolInfo(MolInfo *mi);
int calcSelfSimilarity(MolInfo *mi, Params *parms, int *blev_mat);
float *getDiagonalElements_SS(MolInfo *mi, Params *parms, int *blev_mat); 
int	outputResults(int cid1, int cid2, float ss1, float ss2, float sim, bool S, bool L, bool J, bool M, bool A);
int Simple2Molecules(char *fname1, char *fname2, Params *parms, bool S, bool L, bool J, bool M, bool A, char *iwtsFname);
int MoleculeAgainstList(char *fname1, char *fname2, Params *parms, bool S, bool L, bool J, bool M, bool A, float min_mol_wt, float max_mol_wt, int *filt_db, char *iwtsFname, bool silence);
int ListAgainstItself(char *fname, Params *parms, bool S, bool L, bool J, bool M, bool A, float min_mol_wt, float max_mol_wt, int *filt_in, int *filt_db, bool silence); 
int ListAgainstList(char *fname1, char *fname2, Params *parms, bool S, bool L, bool J, bool M, bool A, float min_mol_wt, float max_mol_wt, int *filt_in, int *filt_db, bool silence); 
float *readAtomWeights(char *fname, int natoms); 



#ifdef __cplusplus
}
#endif



#endif


