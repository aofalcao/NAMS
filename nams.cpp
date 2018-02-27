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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nams.h"
#include "hungarian.h"

#define NUM_ELEMS 11


int MAX_BONDS=0;

SimBond *sbs;			
bool *elims1, *elims2;	
int **the_mat;  //bond similarity matrix
bool BLANK_READ=false;

extern bool USE_MW;
extern float THRES_JACC;


int LU_ELEMS[] ={0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 0, 0, 0, 0, 0, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     0, 8, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

float ADM0[NUM_ELEMS][NUM_ELEMS] = { { 0.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f  } };

float ADM1[11][11] = { { 0.00f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f ,  0.90f  },
 { 0.90f,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.90f ,  0.00f  } };

float ADM2[11][11] = { { 0.00f,  0.48f ,  0.64f ,  0.80f ,  0.96f ,  0.65f ,  0.81f ,  0.97f ,  0.98f ,  0.98f ,  1.00f  },
 { 0.48f,  0.00f ,  0.16f ,  0.50f ,  0.48f ,  0.70f ,  0.50f ,  0.48f ,  0.90f ,  0.50f ,  0.52f  },
 { 0.64f,  0.16f ,  0.00f ,  0.20f ,  0.50f ,  0.40f ,  0.80f ,  0.60f ,  0.70f ,  0.70f ,  0.80f  },
 { 0.80f,  0.50f ,  0.20f ,  0.00f ,  0.16f ,  0.17f ,  0.90f ,  0.17f ,  0.01f ,  0.21f ,  0.27f  },
 { 0.96f,  0.48f ,  0.50f ,  0.16f ,  0.00f ,  0.33f ,  0.17f ,  0.07f ,  0.90f ,  0.14f ,  0.21f  },
 { 0.65f,  0.70f ,  0.40f ,  0.17f ,  0.33f ,  0.00f ,  0.50f ,  0.32f ,  0.07f ,  0.33f ,  0.35f  },
 { 0.81f,  0.50f ,  0.80f ,  0.90f ,  0.17f ,  0.50f ,  0.00f ,  0.16f ,  0.80f ,  0.17f ,  0.21f  },
 { 0.97f,  0.48f ,  0.60f ,  0.17f ,  0.07f ,  0.32f ,  0.16f ,  0.00f ,  0.85f ,  0.07f ,  0.14f  },
 { 0.98f,  0.90f ,  0.70f ,  0.01f ,  0.90f ,  0.07f ,  0.80f ,  0.85f ,  0.00f ,  0.80f ,  0.90f  },
 { 0.98f,  0.50f ,  0.70f ,  0.21f ,  0.14f ,  0.33f ,  0.17f ,  0.07f ,  0.80f ,  0.00f ,  0.07f  },
 { 1.00f,  0.52f ,  0.80f ,  0.27f ,  0.21f ,  0.35f ,  0.21f ,  0.14f ,  0.90f ,  0.07f ,  0.00f  } };

float ADM3[11][11] = { { 0.00f,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  0.00f ,  0.80f ,  0.90f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  0.80f ,  0.00f ,  1.00f ,  1.00f ,  0.70f ,  1.00f ,  1.00f ,  0.90f ,  1.00f ,  1.00f  },
 { 1.00f,  0.90f ,  1.00f ,  0.00f ,  0.80f ,  1.00f ,  0.50f ,  0.90f ,  0.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  0.80f ,  0.00f ,  1.00f ,  1.00f ,  0.10f ,  1.00f ,  0.20f ,  0.30f  },
 { 1.00f,  1.00f ,  0.70f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  0.50f ,  1.00f ,  1.00f ,  0.00f ,  0.90f ,  1.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  0.90f ,  0.10f ,  1.00f ,  0.90f ,  0.00f ,  1.00f ,  0.10f ,  0.20f  },
 { 1.00f,  1.00f ,  0.90f ,  0.00f ,  1.00f ,  1.00f ,  1.00f ,  1.00f ,  0.00f ,  1.00f ,  1.00f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  0.20f ,  1.00f ,  1.00f ,  0.10f ,  1.00f ,  0.00f ,  0.10f  },
 { 1.00f,  1.00f ,  1.00f ,  1.00f ,  0.30f ,  1.00f ,  1.00f ,  0.20f ,  1.00f ,  0.10f ,  0.00f  } };

float ADM4[11][11] = { { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  },
 { 0.00f,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f ,  0.00f  } };

 

float *readAtomWeights(char *fname, int natoms) 
{
	// reads the atomic weights by default for options 0, 1 and 4
	// 
	char seps[]   = " ,\t\n";
	FILE *fil;
	float *wts;
	char line[MAX_LIN_SIZ+1];
	char *tok;
	//puts(fname);
	wts=(float *)malloc(sizeof(float) *natoms);
	if((fil = fopen( fname, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname );
		exit(1);
	}
	fgets(line, MAX_LIN_SIZ, fil);
	fclose(fil);
	//puts(line);
	tok=strtok(line, seps);
	wts[0]=(float)atof(tok);
	//printf("%d --> %f\n", 0, wts[0]);
	for(int a=1; a<natoms; a++) {
		wts[a]=(float)atof(strtok(NULL, seps));
		//printf("%d --> %f\n", a, wts[a]);
	}
	return wts;
}



bool cid_in_flist(int cid, int *filt_list)
{
	//this is the cid filter check. returns true if the cid is in the list and false otherwise
	int count;
	//if the filter is not active, get out!
	if(filt_list==NULL) return true;
	
	count=filt_list[0];
	for(int i=1; i<=count; i++) 
		if(cid==filt_list[i]) return true;
	return false;
}


int setElemsDists(Params *parms, int adm) {
	//Set the atom distance matrix
	parms->LU_ELEMS=&(LU_ELEMS[0]);
	parms->ELEMS_DISTS=(float **)malloc(sizeof(float *)*NUM_ELEMS);
	for(int i=0;i<NUM_ELEMS; i++) {
	    switch(adm) {
	       case 0:
	           parms->ELEMS_DISTS[i]=&(ADM0[i][0]);
	           break;
	       case 1:
	           parms->ELEMS_DISTS[i]=&(ADM1[i][0]);
	           break;
	       case 2:
	           parms->ELEMS_DISTS[i]=&(ADM2[i][0]);
	           break;
	       case 3:
	           parms->ELEMS_DISTS[i]=&(ADM3[i][0]);
	           break;
	       default:
	           parms->ELEMS_DISTS[i]=&(ADM4[i][0]);
	           break;
	    }
	}

	//for(int i=0;i<NUM_ELEMS; i++) {
	//	for(int j=0;j<NUM_ELEMS; j++) {
	//		printf(" %d, %d ->  %f" ,i,j, parms->ELEMS_DISTS[i][j]);
	//	}
	//}
    return 1;
}

int setParams(Params *parms, int ADM, float BS_ALPHA,  float ANRINGS_FAC,  float ACHIR_FAC, 
			  float DBSTEREO_FAC, float BRING_FAC, float BAROM_FAC, float BORDER_FAC, float PEN)
{
  //char fname[1024];
    parms->BS_ALPHA = BS_ALPHA;
    parms->ANRINGS_FAC = ANRINGS_FAC;
    parms->ACHIR_FAC = ACHIR_FAC ;    
    parms->DBSTEREO_FAC = DBSTEREO_FAC;
    parms->BRING_FAC = BRING_FAC;     
    parms->BAROM_FAC = BAROM_FAC;     
    parms->BORDER_FAC = BORDER_FAC;    
    parms->PEN = PEN;
	//printf("ADM: %d", ADM);
	//sprintf(fname, "nelems_dist%01d.txt", ADM);
    setElemsDists(parms, ADM);
    return 1;
}

int readMolInfo( FILE *fil, MolInfo *mi) 
{  
	//if blank_read==true then do not reserve memory
	//this will get a string of data (or a FILE?) and it will parse it writing it to the MolInfo Struct 
	char line[MAX_LIN_SIZ];
	char mname[MAX_LIN_SIZ];
	int natoms, nbonds, nabats;
	char  *ntok=NULL;
	char seps[]   = " ,\t\n";
	int line_aba[ABA_ELEMS];
	int i=0, v;
	

	//0. read the ID
	fgets(line, MAX_LIN_SIZ, fil);
	if(feof(fil)) return feof(fil);
	mi->cid = atoi(strtok(line, seps));
	//this next line will only make sense if the molar weight is included within the file
	if(USE_MW==true) {
		mi->molwt = atoi(strtok(NULL, seps));
	}
	strcpy(mname, strtok(NULL, seps));

	mname[31] = (char)0; //this is very ugly. I just trim the string into the max name size!
	strcpy(mi->mname, mname);
	//1. read the fundamental information
	fgets(line, MAX_LIN_SIZ, fil);
	natoms= atoi(strtok(line, seps));
	nbonds= atoi(strtok(NULL, seps));
	nabats= atoi(strtok(NULL, seps));
	//printf("natoms: %d nbonds: %d naba types: %d\n", natoms, nbonds, nabats);
	mi->natoms=natoms;
	mi->nbonds=nbonds;
	mi->naba_types=nabats;
	if(!BLANK_READ) {
		mi->weights = (float *)malloc(natoms*sizeof(float));
		mi->aba_types = (AbaBond *)malloc(nabats*sizeof(AbaBond));
		mi->abas=(int *)malloc(sizeof(int)*natoms*nbonds);
		mi->levels=(int *)malloc(sizeof(int)*natoms*nbonds);
	}
	//2. read the aba_bond types
	for(int aba=0; aba<nabats; aba++) {
		fgets(line, MAX_LIN_SIZ, fil);
		if(!BLANK_READ) {
			//line_aba[0] = atoi(strtok_s(line, seps, &ntok)); 
			line_aba[0] = atoi(strtok(line, seps)); 
			for(int elem=1; elem< ABA_ELEMS; elem++) {
				//line_aba[elem] = atoi(strtok_s(NULL, seps, &ntok)); 
				line_aba[elem] = atoi(strtok(NULL, seps)); 
			}
			memcpy(&(mi->aba_types[aba]),line_aba, sizeof(int)*ABA_ELEMS);
		}
	}
	
	//3. read the atom x bonds -Bondtypes
	for(int a=0; a<natoms; a++) {
		fgets(line, MAX_LIN_SIZ, fil);
		if(!BLANK_READ) {
			//v= atoi(strtok_s(line, seps, &ntok)); 
			v= atoi(strtok(line, seps)); 
			mi->abas[a*nbonds+0]=v;
			//printf("BondTypes for atom %d:  %d", a, v);
			for(int b=1; b< nbonds; b++) {
				//v= atoi(strtok_s(NULL, seps, &ntok));
				v= atoi(strtok(NULL, seps));
				mi->abas[a*nbonds+b]=v;
				//printf(" %d",mi->abas[a*natoms+b]);
			}
			//printf("\n");
		}
	}
	
	//printf("\n");
	//4. read the atom x bonds -levels
	for(int a=0; a<natoms; a++) {
		fgets(line, MAX_LIN_SIZ, fil);
		if(!BLANK_READ) {
			//v= atoi(strtok_s(line, seps, &ntok)); 
			v= atoi(strtok(line, seps)); 
			mi->levels[a*nbonds+0]=v;
			//printf("Levels for atom %d:  %d", a, mi->levels[a*natoms+0]);
			for(int b=1; b< nbonds; b++) {
				//v= atoi(strtok_s(line, seps, &ntok)); 
				v= atoi(strtok(NULL, seps));
				mi->levels[a*nbonds+b]=v;
				//printf(" %d",mi->levels[a*natoms+b]);
			}
			//printf("\n");
		}
	}
	return 0;
}


int blankReadMolFile( FILE *fil) 
{  
	//this function will not read anything. Its purpose is only for counting molecules and setting max_bonds
	char line[MAX_LIN_SIZ];
	int natoms, nbonds, nabats;
	char  *ntok=NULL;
	char seps[]   = " ,\t\n";
	int i=0, num_mols=0, eof=false;
	//puts("here we are!");
	while(!eof) { 
		//0. read the ID
		fgets(line, MAX_LIN_SIZ, fil);
		//printf(line);
		if(feof(fil)) return num_mols;
		//1. read the fundamental information
		fgets(line, MAX_LIN_SIZ, fil);
		//natoms= atoi(strtok_s(line, seps, &ntok));
		//nbonds= atoi(strtok_s(NULL, seps, &ntok));
		//nabats= atoi(strtok_s(NULL, seps, &ntok));
		natoms= atoi(strtok(line, seps));
		nbonds= atoi(strtok(NULL, seps));
		nabats= atoi(strtok(NULL, seps));
		//printf(" -> %d %d %d\n", natoms, nbonds, nabats);
		if(nbonds>MAX_BONDS) MAX_BONDS=nbonds;
		//2. read the aba_bond types
		for(int aba=0; aba<nabats; aba++) fgets(line, MAX_LIN_SIZ, fil);	
		//3. read the atom x bonds -Bondtypes
		for(int a=0; a<natoms; a++) fgets(line, MAX_LIN_SIZ, fil);
		//4. read the atom x bonds -levels
		for(int a=0; a<natoms; a++) fgets(line, MAX_LIN_SIZ, fil);
		num_mols++;
		//printf("NMols: %5d  natoms: %d nbonds: %d naba types: %d\n", num_mols, natoms, nbonds, nabats); 
		eof=feof(fil);
	}
	return num_mols;
}




int showMolInfo(MolInfo *mi) {
	//this is a check function to verify we manage to read a given molecule correctly
	int natoms, nbonds, nabats, v; 
	natoms= mi->natoms;
	nbonds=mi->nbonds;
	nabats=mi->naba_types;
	int myaba[ABA_ELEMS];
	printf("natoms: %d nbonds: %d naba types: %d\n", natoms, nbonds, nabats);

	//1. show the aba_bond types
	for(int aba=0; aba<nabats; aba++) {
		printf("ABA_Bond Info %d: ", aba);
		memcpy(myaba, &(mi->aba_types[aba]), sizeof(int)*ABA_ELEMS);
		for(int elem=0; elem< ABA_ELEMS; elem++) {
			printf(" %d",myaba[elem]);
		}
		printf("\n");
	}
	printf("\n");
	//2. read the atom x bonds -Bondtypes
	for(int a=0; a<natoms; a++) {
		printf("BondTypes for atom %d: ", a);
		for(int b=0; b< nbonds; b++) {
			v=mi->abas[a*nbonds+b];
			printf(" %d",v);
		}
		printf("\n");
	}
	printf("\n");
	//4. read the atom x bonds -levels
	for(int a=0; a<natoms; a++) {
		printf("Levels for atom %d: ", a);
		for(int b=0; b< nbonds; b++) {
			v=mi->levels[a*nbonds+b];
			printf(" %d",v);
		}
		printf("\n");
	}
	return 0;
}

int readSimAtomMAtrix(char *fname) 
{
	//here we read an atom similarity matrix, and will put it in a integer matrix
	return 1;
}


int compare_aba_bonds(AbaBond *aba1, AbaBond *aba2, Params *parms) 
{

	float simil=1.0;
	//int numb1, nrings1, chir1, numb2, nrings2, chir2, inring, arom, order, dbcistrans;
	int a11, a12, a21, a22;
	//lookup the atom indexes for checking similarity
	a11=parms->LU_ELEMS[aba1->numb1];
	a12=parms->LU_ELEMS[aba1->numb2];
	a21=parms->LU_ELEMS[aba2->numb1];
	a22=parms->LU_ELEMS[aba2->numb2];
	//now we can compute the atom similarity
    simil *= (1.0f - parms->ELEMS_DISTS[a11][a21]);
	//printf("%d %d %f simil: %f\n", a11, a21, parms->ELEMS_DISTS[a11][a21], simil);
    simil *= (1.0f - parms->ELEMS_DISTS[a12][a22]);
	//printf("%d %d %f simil: %f\n", a12, a22, parms->ELEMS_DISTS[a12][a22], simil);

    if(aba1->chir1*aba2->chir1 ==-1) simil *= parms->ACHIR_FAC;
	//printf("CHIR1 %d %d simil: %f\n", aba1->chir1, aba2->chir1, simil);
    if (aba1->chir2*aba2->chir2 ==-1) simil *= parms->ACHIR_FAC;
	//printf("CHIR2 %d %d simil: %f\n", aba1->chir2, aba2->chir2, simil);

	//simil *= (1.0f - abs(aba1->nrings1-aba2->nrings1)*parms->ANRINGS_FAC/4.0f);
    //simil *= (1.0f - abs(aba1->nrings2-aba2->nrings2)*parms->ANRINGS_FAC/4.0f);
	//this is a really important correction! (11/20/2017)
	simil *= powf(parms->ANRINGS_FAC, (float)abs(aba1->nrings1-aba2->nrings1));
	simil *= powf(parms->ANRINGS_FAC, (float)abs(aba1->nrings2-aba2->nrings2));

    if(aba1->inring != aba2->inring) simil *= parms->BRING_FAC;
    if(aba1->arom  != aba2->arom) simil   *= parms->BAROM_FAC;
    if(aba1->order != aba2->order) simil  *= parms->BORDER_FAC;
    if(aba1->dbcistrans * aba2->dbcistrans == -1) simil *= parms->DBSTEREO_FAC;
        
    //return (int)((2.0f*simil-1.0f)*1000);
	return (int)(simil*100);
}

int *getAbaBondsCompMatrix(MolInfo *mi1, MolInfo *mi2, Params *parms) 
{
	//this function receives two molecules and will produce a matrix with the comparison of all bonds one agains the other
	//using the order of the ababonds in each atom 
	// the matrix is represented as an array of integers(?) to access each the formula is
	//mat[a2+na2*a1]
	int *mat=NULL;
	int res;

	mat=(int *)malloc(mi1->naba_types*mi2->naba_types*sizeof(int));
	//memset(mat,1, mi1->naba_types*mi2->naba_types*sizeof(int)); //LIXO!!!!!!!!!!
	
	//res=compare_aba_bonds(&(mi1->aba_types[2]), &(mi1->aba_types[4]), parms); 
	//printf("%6.3f ", res);
	for(int i=0; i<mi1->naba_types; i++) {
		for(int j=0; j<mi2->naba_types; j++) {
			res=compare_aba_bonds(&(mi1->aba_types[i]), &(mi2->aba_types[j]), parms); 
			mat[mi2->naba_types*i+ j]=res;
			//printf("%6d ", res);
		}
		//puts("");
	}
	
	return mat;
}


int *calcBondLevelsMatrix(Params *parms) 
{
	//this is a very important function and for efficiency should be computed only ONCE and the results stored in one array
	//the idea is that we should compute the floating point stuf once and only once
	// to make everything really fast, the striangular matrix will be symmetrical
	int *mat=(int *)malloc(sizeof(int)*MAX_LEVELS*MAX_LEVELS);
	int v, mx;
	for(int i=0; i<MAX_LEVELS-1; i++) {
		for(int j=i; j<MAX_LEVELS; j++) {
			//v=(int)(100*1.0f/powf((abs(i-j)+j+1.0f), parms->BS_ALPHA));
			//modified at 20/11/2017
			mx=i;
			if(parms->BS_ALPHA>0.0f) {
				if(j>mx) mx=j;
				v=(int)(100.0*powf(parms->BS_ALPHA, (float)(mx+abs(i-j))));
			}else {
				if(i==j) v=100; else v=0;
			}
			mat[i*MAX_LEVELS+j] = v;
			mat[j*MAX_LEVELS+i] = v;
		}
	}
	return mat;
}





int simbond_compare( const void *sim1, const void *sim2 )
{
	// this is a callback for qsort, for sorting bond similaritis in SimBonds;
	int s1= ((SimBond *)sim1)->sim;
	int s2= ((SimBond *)sim2)->sim;
	if(s1 > s2) return -1;
	if(s1 < s2) return 1;
	return 0;
}

 


int matchBonds(int a1, int a2, MolInfo *mi1, MolInfo *mi2, int *aba_ms, int *blev_mat) 
{
	// this is one of the MOST CRITICAL functions of this whole story. Must be carefully crafted
	// here we will make a matrix with as many rows as the number of bonds of mol1 and rows as n of bonds for mol2
	// the matrix will have the comparison score USING the levels for each bond of mol1 against each bond of mol2
	// aba_ms -> similarity matrix between aba_bond types
	// blev_mat -> levels matrix - >precomputed factors that account for the level distances between bonds
	// a1 and a2 -> atom indexes for getting the corresponding rows out of the molInfos

	int *levs1, *levs2, *abas1, *abas2;
	int lev1, lev2, abat1, abat2, nfac = mi2->naba_types;
	int sim=0, mscore=0, ce1=0, ce2=0;
	//thes variables are just for information purposes of counting the best possible matches on rows and columns
	int sbr=0, sbc=0, max_sim=0;
	int nbonds1, nbonds2, max_row=0, max_col=0;
	//float minrc, maxrc, arow, acol, sum_sims=0, sum_sims2=0, avg_sims=0, var_sims=0, simf=0;
	
	//SimBond *sbs;  //this probably could be allocated globally as for the levels matrix and used over and over again for all molecules
	//bool *elims1, *elims2; //for these two the same caveat applies;
	nbonds1=mi1->nbonds;
	nbonds2=mi2->nbonds;
	if(MAX_BONDS==0) {
		elims1=(bool *)malloc(sizeof(bool)* nbonds1);
		elims2=(bool *)malloc(sizeof(bool)* nbonds2);
		the_mat =(int **)malloc(sizeof(int *)*nbonds1);
		for(int i=0; i<nbonds1; i++) the_mat[i]=(int *)malloc(sizeof(int)*nbonds2);
		//the_mat =(int *)malloc(sizeof(int *)*nbonds1*nbonds2);

		//reset elims
		memset(elims1, false, sizeof(bool)*nbonds1);
		memset(elims2, false, sizeof(bool)*nbonds2);
	} else {
		memset(elims1, false, sizeof(bool)*MAX_BONDS);
		memset(elims2, false, sizeof(bool)*MAX_BONDS);
	}

	abas1 = &(mi1->abas[nbonds1*a1]);
	abas2 = &(mi2->abas[nbonds2*a2]);
	levs1 = &(mi1->levels[nbonds1*a1]);
	levs2 = &(mi2->levels[nbonds2*a2]);


	//this is it! Compute the bond matrix against 2 atoms
	//printf("bond matrix against 2 atoms %d %d\n",mi1->nbonds, mi2->nbonds );
	for(int i=0; i<nbonds1; i++) {
		abat1 = abas1[i];
		lev1  = levs1[i];
		//printf("----->");
		for(int j=0; j<nbonds2; j++) {
			abat2 = abas2[j];
			lev2  = levs2[j];
			sim=blev_mat[lev1*MAX_LEVELS+lev2] * aba_ms[abat1*nfac+abat2];
			if(sim>max_sim) {
				max_sim = sim;
				max_row=i;
				max_col=j;
			}
			//the_mat[i*nbonds2+j] = sim;
			the_mat[i][j]=sim;
		}
	}


	

	elims1[max_row]=true;
	elims2[max_col]=true;
	mscore=max_sim;
	//memset(the_mat[max_row], -1, sizeof(int)*nbonds2);
	//for(int i=0; i<nbonds1; i++) the_mat[i][max_col]=-1;

	//this could speed up by counting the number of hits and breaking when hits = nbonds 
	while(ce1<nbonds1 && ce2 <nbonds2) {
		max_sim=-1;
		for(int row=0; row< nbonds1; row++) {
			if(!elims1[row]) {
				for(int col=0; col<nbonds2; col++) {
					if(!elims2[col] && the_mat[row][col] > max_sim) {
						//max_sim=the_mat[row][col];
						max_sim=the_mat[row][col];
						max_row=row;
						max_col=col;
					}
				}
			}
		}		
		elims1[max_row]=true;
		elims2[max_col]=true;
		mscore+=max_sim;
		//memset(the_mat[max_row], -1, sizeof(int)*nbonds2);
		//for(int i=0; i<nbonds1; i++) the_mat[i][max_col]=-1;
		ce1++;
		ce2++;
	}



	//-----------------------------END TODO
	if(MAX_BONDS==0) {
		free(sbs);
		free(elims2);
		free(elims1);
		for(int i=0;i<nbonds1; i++) free(the_mat[i]);
		free(the_mat);
	}
	return mscore;
}



int *calcAtomMatchingMatrix(MolInfo *mi1, MolInfo *mi2, int *aba_match_scores, int *blev_mat, float pen)
{
	//this is a critical function. It will compare all atoms on mol1 to all atoms of mol2 and create a score for each match
	//mat will hold the scores
	int*  mat=(int*)malloc(mi1->natoms*mi2->natoms*sizeof(int));
	int c=0, bscore;
	for(int i=0; i<mi1->natoms;i++) {
		//VERY VERY IMPORTANT! WHEN WEIGTHING THE ATOMS IMPORTANCE THIS IS WHERE WE SHOULD ACCOUNT FOR THEM
		for(int j=0; j < mi2->natoms; j++) {
			bscore = matchBonds(i, j, mi1, mi2, aba_match_scores, blev_mat); 
			//the penalties stuff
			bscore -= (int)(pen*abs(mi1->nbonds-mi2->nbonds)*100.0);
			mat[c]=bscore;
			c++;
			//printf("%5.2f ", bscore/10000.0f);
		}
		//puts("");
	}
	return mat;
}


int nams_runner(MolInfo *mi1, MolInfo *mi2, Params *parms, int *blev_mat, bool M, bool A, float *wts)
{
	//Here is the backbone of it all. 2 molecules to be compared 
	//d) the atom matching matrix (M)
	//e) the alignment produced (A)

	int	*aba_match_scores; //the scores for each aba_bond match between 2 molecules
	int	*atom_matrix;
	int	col;
	int	v, final_score=0;
	hungarian_problem_t p;

	aba_match_scores = getAbaBondsCompMatrix(mi1, mi2, parms); 
	
	//for(int i=0; i<mi->naba_types; i++) {
	//	for(int j=0; j<mi->naba_types; j++) printf("%6d ", aba_match_scores[mi->naba_types*i+ j]);
	//	puts("");
	//}	

	atom_matrix=calcAtomMatchingMatrix(mi1, mi2, aba_match_scores, blev_mat, parms->PEN);
	//puts("");

	
	int** m = array_to_matrix(atom_matrix, mi1->natoms, mi2->natoms); //THIS is where the weights will strike in!
	if(wts) {
		for(int row=0; row<mi1->natoms; row++) {
			for(int col=0; col<mi2->natoms; col++) {
				if(m[row][col]>0) m[row][col]=(int)(powf(m[row][col]/10000.0f, wts[row])*10000);
				else m[row][col] = 0;
			}
		}
	}
	
	int matrix_size = hungarian_init(&p, m , mi1->natoms, mi2->natoms, HUNGARIAN_MODE_MAXIMIZE_UTIL) ;

	//here it is!!!
	hungarian_solve(&p);
	//this is just for writing
	if(M) {
		puts("Atom matching matrix");
		printf("    ");
		for(int j =0; j<mi2->natoms; j++) printf("%6d ", j+1); 
		puts("");
		for(int i =0; i<mi1->natoms; i++) {
			printf("%3d ", i+1);
			for(int j =0; j<mi2->natoms; j++) {
				printf("%6.2f ", atom_matrix[i*mi2->natoms+j]/10000.0f); 
			}
			puts("");
		}
	}
	//this is where the actual result is calculated
	if(A) puts("Atom matching and scores");
	for (int row = 0; row < mi1->natoms; row++) {
		col=p.my_assig[row];
		if(col<mi2->natoms && row < mi1->natoms) {
			//v=atom_matrix[row*mi2->natoms+col]; // /10000.0f;
			v=m[row][col];
			if(A) printf("%3d %3d -> %5.2f\n", row+1, col+1,  v/10000.0f);
			final_score+=v;
			//printf("\t%d %d -> %d\n", row, col , v);
		}
	}
	//printf("Matching Score: %7.3f\n", final_score);
	hungarian_free(&p);
	free_array_to_matrix(m, mi1->natoms); 
	free(atom_matrix);
	free(aba_match_scores);
	return final_score;
}

int freeMolInfo(MolInfo *mi) 
{
	free(mi->weights);
	free(mi->aba_types);;
	free(mi->abas);
	free(mi->levels);
	return 1;
}


int calcSelfSimilarity(MolInfo *mi, Params *parms, int *blev_mat) 
{
	//this function returns the self similarity in one molecule by simply compu7ting the diagonal elements
	int*  mat=(int*)malloc(mi->natoms*mi->natoms*sizeof(int));
	int bscore=0;
	int  *aba_match_scores;

	aba_match_scores = getAbaBondsCompMatrix(mi, mi, parms); 

	for(int i=0; i<mi->natoms;i++) {
		//printf("ATOM: %d\n", i);
		bscore += matchBonds(i, i, mi, mi, aba_match_scores, blev_mat); 
		//printf("%d %d\n", i, bscore);
	}
	free(mat);
	free(aba_match_scores);
	return bscore;
}

float *getDiagonalElements_SS(MolInfo *mi, Params *parms, int *blev_mat) 
{
	//this function returns the diagonal elements of the atom matrix from self similary comparison
	int*	mat=(int*)malloc(mi->natoms*mi->natoms*sizeof(int));
	float   *diag=(float *)malloc(mi->natoms*sizeof(float));
	int		*aba_match_scores;

	aba_match_scores = getAbaBondsCompMatrix(mi, mi, parms); 

	for(int i=0; i<mi->natoms;i++) {
		diag[i] =matchBonds(i, i, mi, mi, aba_match_scores, blev_mat)/10000.0f; 
		//printf("\t\t%d %7.4f\n", i, diag[i]);
	}
	free(mat);
	free(aba_match_scores);
	return diag;
}



int	outputResults(int cid1, int cid2, float ss1, float ss2, float sim, bool S, bool L, bool J, bool M, bool A)
{
	float jacc=sim/(ss1+ss2-sim);
	if(A || M) {
		if(L || S || J) puts("");
		if(L)  {
			printf("Self Similarity: %4d -> %7.3f\n", cid1, ss1);
			printf("Self Similarity: %4d -> %7.3f\n", cid2, ss2);
		}
		if(S) printf("Similarity: (%d %d) -> %7.3f\n", cid1, cid2, sim);
		if(J) printf("Jaccard Score: (%d %d) -> %7.5f\n", cid1, cid2, sim/(ss1+ss2-sim));
	} else {
		if(jacc>=THRES_JACC) {
			if(L || S || J) {
				printf("%10d%10d%", cid1, cid2);
				if(L) printf(" %7.3f %7.3f", ss1, ss2); else printf("                ");
				if(S) printf(" %7.3f", sim); else printf("        ");
				if(J) printf(" %7.4f",sim/(ss1+ss2-sim)); else printf("        ");
				puts("");
			}
		}
	}
	fflush(stdout);
	return 1;
}



int Simple2Molecules(char *fname1, char *fname2, Params *parms, bool S, bool L, bool J, bool M, bool A, char *iwtsFname) 
{
/*
a) their similarity scores (S)
b) their self similarity scores (L)
c) their Jaccard coefficient (J)
d) the atom matching matrix (M)
e) the alignment produced (A)
*/
	MolInfo		*mi1, *mi2;
	FILE		*fil1, *fil2;
	int			*blev_mat;			//the precomputed level factors
	float		ss1=0, ss2, sim;
	float		*wts;
	float		*ss1_diags;

	blev_mat = calcBondLevelsMatrix(parms);

	mi1=(MolInfo *)malloc(sizeof(MolInfo));
	mi2=(MolInfo *)malloc(sizeof(MolInfo));
	
	if((fil1 = fopen( fname1, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname1 );
		exit(1);
	}
	readMolInfo(fil1, mi1);
	fclose(fil1);
	if((fil2 = fopen( fname2, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname2 );
		exit(1);
	}
	if(iwtsFname[0]){
		wts = readAtomWeights(iwtsFname, mi1->natoms); 
		ss1_diags = getDiagonalElements_SS(mi1, parms, blev_mat);
		for(int r=0; r<mi1->natoms; r++) ss1 +=(float)powf(ss1_diags[r], wts[r]);
	} else {
		ss1 = calcSelfSimilarity(mi1, parms, blev_mat)/10000.0f;
	}



	readMolInfo(fil2, mi2);
	fclose(fil2);

	ss2 = calcSelfSimilarity(mi2, parms, blev_mat)/10000.0f;
	sim = nams_runner(mi1, mi2, parms, blev_mat, M, A, wts)/10000.0f;
	puts("");
	if(L) {
		printf("Self Similarity: %4d -> %7.3f\n", mi1->cid, ss1);
		printf("Self Similarity: %4d -> %7.3f\n", mi2->cid, ss2);
	}
	
	if(S) printf("Similarity: (%d %d) -> %7.3f\n", mi1->cid, mi2->cid, sim);
	if(J) printf("Jaccard Score: (%d %d) -> %7.5f\n", mi1->cid, mi2->cid, sim/(ss1+ss2-sim));
	freeMolInfo(mi1); 
	freeMolInfo(mi2); 
	free(mi1);
	free(mi2);

	free(blev_mat);
	return 1;
}


int MoleculeAgainstList(char *fname1, char *fname2, Params *parms, bool S, bool L, bool J, bool M, bool A, 
						float min_mol_wt, float max_mol_wt, int *filt_db, char *iwtsFname, bool silence) 
{
/*
a) their similarity scores (S)
b) their self similarity scores (L)
c) their Jaccard coefficient (J)
d) the atom matching matrix (M)
e) the alignment produced (A)
*/
	MolInfo *mi1, *mi2;
	FILE    *fil1, *fil2;
	int     *blev_mat;			//the precomputed level factors
	float	ss1=0, ss2, sim;
	bool	goforit;
	int		eof=false;
	float	*ss1_diags;
	float	*wts=NULL;
	int		ncomp=0;
	
	//puts("step1");
	blev_mat = calcBondLevelsMatrix(parms);
	 
	mi1=(MolInfo *)malloc(sizeof(MolInfo));
	mi2=(MolInfo *)malloc(sizeof(MolInfo));


	//puts("step2");	
	if((fil1 = fopen( fname1, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname1 );
		exit(1);
	}
	eof=readMolInfo(fil1, mi1);
	fclose(fil1);
	//printf("Weight of mol1: %d\n", mi1->molwt);
	//The base molecule has been read. Now, if there are input weights, go for it 
	if(iwtsFname[0]){
		wts = readAtomWeights(iwtsFname, mi1->natoms); 
		ss1_diags = getDiagonalElements_SS(mi1, parms, blev_mat);
		for(int r=0; r<mi1->natoms; r++) ss1 +=(float)powf(ss1_diags[r], wts[r]);
	} else {
		ss1 = calcSelfSimilarity(mi1, parms, blev_mat)/10000.0f;
	}

	//puts("step4");	
	//fil2=fopen(fname2, "rt");
	if((fil2 = fopen( fname2, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname2 );
		exit(1);
	}
	//here it starts
	eof=readMolInfo(fil2, mi2);
	//printf("Weight of mol2: %d\n", mi2->molwt);
	//puts("step5");	
	if((L || S || J) && (!A && !M)) printf("    molid1      molid2     ss1     ss2   sim12   Sscore\n");
	while(!eof) {
		if(USE_MW==true) {
			//printf("--> %d   -   %f %f\n", mi2->molwt , max_mol_wt*mi1->molwt, min_mol_wt*mi1->molwt);
			if(mi2->molwt <= max_mol_wt*mi1->molwt && mi2->molwt >= min_mol_wt*mi1->molwt) {
				goforit=true;
				if(filt_db) goforit=cid_in_flist(mi2->cid, filt_db);
				if(goforit==true) {
					if(L || J) ss2 = calcSelfSimilarity(mi2, parms, blev_mat)/10000.0f;
					sim = nams_runner(mi1, mi2, parms, blev_mat, M, A, wts)/10000.0f;
					outputResults(mi1->cid, mi2->cid, ss1, ss2, sim, S, L, J, M, A);
					freeMolInfo(mi2);
					ncomp++;
				}
			}
		} else {
			goforit=true;
			if(filt_db) goforit=cid_in_flist(mi2->cid, filt_db);
			if(goforit==true) {
				if(L || J) ss2 = calcSelfSimilarity(mi2, parms, blev_mat)/10000.0f;
				sim = nams_runner(mi1, mi2, parms, blev_mat, M, A, wts)/10000.0f;
				outputResults(mi1->cid, mi2->cid, ss1, ss2, sim, S, L, J, M, A);
				ncomp++;
				freeMolInfo(mi2);
			}
		}
		eof=readMolInfo(fil2, mi2);
	}
	if(!silence) fprintf(stderr, "Number of molecular alignments performed: %d\n", ncomp);

	fclose(fil2);
	freeMolInfo(mi1); 
	free(mi1);
	free(mi2);

	free(blev_mat);
	return 1;
	
}


int ListAgainstItself(char *fname, Params *parms, bool S, bool L, bool J, bool M, bool A, float min_mol_wt, float max_mol_wt, int *filt_in, int *filt_db, bool silence) 
{
/*
	a) their similarity scores (S)
	b) their self similarity scores (L)
	c) their Jaccard coefficient (J)
	d) the atom matching matrix (M)
	e) the alignment produced (A)

*/ 
	MolInfo *miList, *mi1, *mi2;
	FILE	*fil;
	int		num_mols;
	float	*selfies;  //self similarities array
	int		*blev_mat;
	float	sim;
	bool	goforit;
	float	*wts=NULL;  //wts should never be used in this function so this NULL for nams_runer
	int		ncomp=0;	//number of comparisons made
	int		start_loop=0, end_loop=0;
	blev_mat = calcBondLevelsMatrix(parms);

	//first get the number of molecules and initialize an internal array
	//fil=fopen(fname, "rt");
	if((fil = fopen( fname, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname );
		exit(1);
	}
	num_mols=blankReadMolFile(fil);
	fclose(fil);
	if(!silence) {
		fprintf(stderr, "Number of molecules in file: %d\n", num_mols);
		fprintf(stderr, "MAX BONDS: %d\n", MAX_BONDS);
	}
	miList=(MolInfo *)malloc(sizeof(MolInfo)*num_mols);

	//initialize the bond arrays to a a max size, so as to avoid heap frafmentation!
	elims1=(bool *)malloc(sizeof(bool)* MAX_BONDS);
	elims2=(bool *)malloc(sizeof(bool)* MAX_BONDS);
	sbs = (SimBond *)malloc(sizeof(SimBond)* MAX_BONDS*MAX_BONDS);
	//the_mat = (int *)malloc(sizeof(int *)*MAX_BONDS*MAX_BONDS);
	the_mat =(int **)malloc(sizeof(int *)*MAX_BONDS);
	for(int i=0; i<MAX_BONDS; i++) the_mat[i]=(int *)malloc(sizeof(int)*MAX_BONDS);



	//reopen the file and actually read the molecules storing them in the miList array
	//fil=fopen(fname, "rt");
	if((fil = fopen( fname, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname );
		exit(1);
	}
	for(int mol=0; mol<num_mols; mol++) {
		readMolInfo(fil, &(miList[mol]));
		//showMolInfo(&(miList[mol])); // !!!!!
	}
	fclose(fil);
	//puts("mols read and done!");

	//now start the calculations. First, compute the self similarities if needed
	if(L || J) {
		selfies=(float *)malloc(sizeof(float)*num_mols); 
		for(int mol=0; mol<num_mols; mol++) {
			selfies[mol]=calcSelfSimilarity(&(miList[mol]), parms, blev_mat)/10000.0f;
		}
	}
	//now iterate through all the molecules
	if((L || S || J) && (!A && !M)) printf("    molid1      molid2     ss1     ss2   sim12   Sscore\n");
	if(filt_in) end_loop=0; else end_loop=1;
	for(int mol1=0; mol1<num_mols-end_loop; mol1++) {
		if(filt_in) start_loop=0; else start_loop=mol1+1;
		//printf("%d is the startloop\n", start_loop);
		for(int mol2=start_loop; mol2<num_mols; mol2++) {
			mi1 = &(miList[mol1]);
			mi2 = &(miList[mol2]);
			if(cid_in_flist(mi1->cid, filt_in)) {
				//printf("%d is in fin\n", mi1->cid);
				if(USE_MW==true && mi2->molwt <= max_mol_wt*mi1->molwt && mi2->molwt >= min_mol_wt*mi1->molwt) {
					goforit=true;
					if(filt_db) goforit=cid_in_flist(mi2->cid, filt_db);
					if(goforit==true) {
						//printf("\t%d is in fdb\n", mi2->cid);
						sim = nams_runner(mi1, mi2, parms, blev_mat, M, A, wts)/10000.0f;
						outputResults(mi1->cid, mi2->cid, selfies[mol1], selfies[mol2], sim, S, L, J, M, A);
						ncomp++;
					}
				}
			} //else printf("%d is NOT in fin\n", mi1->cid);
		}
	}
	if(!silence) fprintf(stderr, "Number of molecular alignments performed: %d\n", ncomp);
	//clean up
	free(elims2);
	free(elims1);
	free(sbs);
	for(int i=0;i<MAX_BONDS; i++) free(the_mat[i]);
	free(the_mat);


	return 1;
}


int ListAgainstList(char *fname1, char *fname2, Params *parms, bool S, bool L, bool J, bool M, bool A, float min_mol_wt, float max_mol_wt, int *filt_in, int *filt_db, bool silence)
{
/*
	a) their similarity scores (S)
	b) their self similarity scores (L)
	c) their Jaccard coefficient (J)
	d) the atom matching matrix (M)
	e) the alignment produced (A)

*/
	MolInfo		*miList1, *miList2, *mi1, *mi2;
	FILE		*fil;
	int			num_mols1, num_mols2;
	float		*selfies1, *selfies2;  //self similarities array
	int			*blev_mat;
	float		sim;
	float		*wts=NULL; //weights should NEVER be used in this function, so this is NULL for nams_runner
	int			ncomp=0;
	bool		goforit;
	blev_mat = calcBondLevelsMatrix(parms);

	//first get the number of molecules and initialize an internal array
	//fil=fopen(fname1, "rt");
	if((fil = fopen( fname1, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname1 );
		exit(1);
	}
	num_mols1=blankReadMolFile(fil);
	fclose(fil);
	if(!silence) {
		fprintf(stderr, "Number of molecules in file in: %d\n", num_mols1);
		fprintf(stderr, "MAX BONDS: %d\n", MAX_BONDS);
	}
	miList1=(MolInfo *)malloc(sizeof(MolInfo)*num_mols1);

	//first get the number of molecules and initialize an internal array
	//fil=fopen(fname2, "rt");
	if((fil = fopen( fname2, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname2 );
		exit(1);

	}
	num_mols2=blankReadMolFile(fil);
	fclose(fil);
	if(!silence) {
		fprintf(stderr, "Number of molecules in file db: %d\n", num_mols2);
		fprintf(stderr, "MAX BONDS: %d\n", MAX_BONDS);
	}
	miList2=(MolInfo *)malloc(sizeof(MolInfo)*num_mols2);



	//initialize the bond arrays to a a max size, so as to avoid heap fragmentation!
	elims1=(bool *)malloc(sizeof(bool)* MAX_BONDS);
	elims2=(bool *)malloc(sizeof(bool)* MAX_BONDS);
	sbs = (SimBond *)malloc(sizeof(SimBond)* MAX_BONDS*MAX_BONDS);
	the_mat =(int **)malloc(sizeof(int *)*MAX_BONDS);
	for(int i=0; i<MAX_BONDS; i++) the_mat[i]=(int *)malloc(sizeof(int)*MAX_BONDS);



	//reopen the files and actually read the molecules storing them in the miList arrays
	//fil=fopen(fname1, "rt");
	if((fil = fopen( fname1, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname1 );
		exit(1);
	}
	for(int mol=0; mol<num_mols1; mol++) readMolInfo(fil, &(miList1[mol]));
	fclose(fil);

	//rtheother file
	//fil=fopen(fname2, "rt");
	if((fil = fopen( fname2, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname2 );
		exit(1);
	}
	for(int mol=0; mol<num_mols2; mol++) readMolInfo(fil, &(miList2[mol]));
	fclose(fil);
	//puts("mols read and done!");

	//now start the calculations. First, compute the self similarities if needed
	if(L || J) {
		selfies1=(float *)malloc(sizeof(float)*num_mols1); 
		for(int mol=0; mol<num_mols1; mol++) selfies1[mol]=calcSelfSimilarity(&(miList1[mol]), parms, blev_mat)/10000.0f;
		
		selfies2=(float *)malloc(sizeof(float)*num_mols2); 
		for(int mol=0; mol<num_mols2; mol++) selfies2[mol]=calcSelfSimilarity(&(miList2[mol]), parms, blev_mat)/10000.0f;
	}
	//now iterate through all the molecules
	if((L || S || J) && (!A && !M)) printf("    molid1      molid2     ss1     ss2   sim12   Sscore\n");
	for(int mol1=0; mol1<num_mols1; mol1++) {
		for(int mol2=0; mol2<num_mols2; mol2++) {
			mi1 = &(miList1[mol1]);
			mi2 = &(miList2[mol2]);
			if(cid_in_flist(mi1->cid, filt_in)) {
				if(USE_MW==true && mi2->molwt <= max_mol_wt*mi1->molwt && mi2->molwt >= min_mol_wt*mi1->molwt) {
					goforit=true;
					if(filt_db) goforit=cid_in_flist(mi2->cid, filt_db);
					if(goforit==true) {
						sim = nams_runner(mi1, mi2, parms, blev_mat, M, A, wts)/10000.0f;
						outputResults(mi1->cid, mi2->cid, selfies1[mol1], selfies2[mol2], sim, S, L, J, M, A);
						ncomp++;
					}
				}
			}
		}
	}
	if(!silence) fprintf(stderr, "Number of molecular alignments performed: %d\n", ncomp);
	//clean up
	free(elims2);
	free(elims1);
	free(sbs);
	for(int i=0;i<MAX_BONDS; i++) free(the_mat[i]);
	free(the_mat);
	free(miList1);
	free(miList2);
	free(selfies1);
	free(selfies2);
	return 1;
}

