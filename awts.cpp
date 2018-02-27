#include "awts.h"
#include "nams.h"

extern int MAX_BONDS;	
extern SimBond *sbs;			
extern int **the_mat;  

extern bool *elims1;
extern bool *elims2;	



int WeightsEvaluator(float ***matrices, int **mat, float *wts, int nmols, int natoms1, int *natoms2, float *ss1, float*ss2, float *scores, bool J)
{
	int		matrix_size;
	float	fs_f=0;		//final score as a float
	int		col, row, final_score=0, v;
	float	real_ss1=0.0f; //the real self similarity of molecule 1 when applied the weights
	hungarian_problem_t p;

	//now this is the tricky part
	//now for all molecules compute the matrix values 
	for(int i=0; i<nmols; i++) {
		for(int r=0; r<natoms1; r++) {
			for(int c=0; c<natoms2[i]; c++) {
				//now apply the weights as exponents for each element of each row
				if(matrices[i][r][c]>0) mat[r][c]=(int)(powf(matrices[i][r][c]/10000.0f, wts[r])*10000);
				else mat[r][c] = 0;
				//printf("\t%d %d -> %f %f  -->  %d\n", r, c, wts[r], matrices[i][r][c],  mat[r][c]);
			}
		}
		//now solve the system
		//first create the structure
		matrix_size = hungarian_init(&p, mat, natoms1, natoms2[i], HUNGARIAN_MODE_MAXIMIZE_UTIL) ;
		//here it is!!!
		hungarian_solve(&p);
		final_score=0;
		for (row = 0; row < natoms1; row++) {
			col=p.my_assig[row];
			if(col<natoms2[i] && row < natoms1) {
				v=mat[row][col]; 
				final_score+=v;
				//printf("\t%d %d -> %d\n", row, col , v);
			}
		}
		hungarian_free(&p);
		real_ss1=0.0;
		//the ss1 are the diagonal element, so we must compute the real score as the sum of the POWEREDdiagonal elements 
		for(int r=0; r<natoms1; r++){
			real_ss1+=(float)powf(ss1[r], wts[r]);
			//printf("\t\t%d %7.4f\n", r, ss1[r]);
		}
		fs_f = final_score/10000.0f;
		if(J) scores[i] =fs_f/(real_ss1+ss2[i]-fs_f);
		else scores[i] =fs_f;
		//printf("\t%7.4f %7.4f %7.4f --> %7.4f      %d  %d\n", real_ss1, ss2[i], fs_f, scores[i], natoms1, natoms2[i]);
	}
	
	return 0;
}


int GenerateAtomWeights(char *fname1, char *fname2, Params *parms, bool J, 
						int *filt_db, int nruns, char *wtsFname, char *iwtsFname, float spread) 
{
/*
this function receives a molecule, a database, and a  LIST of molecule IDs present in the database
firstly it will compute the atom matrices and 
*/
	int	*aba_match_scores; //the scores for each aba_bond match between 2 molecules
	int *atom_matrix;
	MolInfo *mi1, *mi2;
	FILE    *fil1, *fil2;
	FILE	*wFile=NULL; 
	int     *blev_mat;			//the precomputed level factors
	float	*ss1, *ss2;			//the self similarities
	bool	goforit;
	int		eof=false;
	int		m;
	int		num_mols;
	float	***matrices;
	int		natoms1, *natoms2, max_natoms2=0;
	float	*scores;
	float	*wts;
	float	*init_wts;
	float	rnd;
	float	sumw=0, fac=0;

	if(wtsFname[0])	wFile=fopen(wtsFname, "wt");
	//get the max number of mols to compare and use the value to allocate the arrays of precomputed data
	if((fil2 = fopen( fname2, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname2 );
		exit(1);
	}
	num_mols=blankReadMolFile(fil2);
	fclose(fil2);
	if(filt_db) num_mols=filt_db[0]; //the  first element is the number of elements on the list

	//initialize the bond arrays to a a max size and keep it to the end, so as to avoid heap frafmentation!
	elims1=(bool *)malloc(sizeof(bool)* MAX_BONDS);
	elims2=(bool *)malloc(sizeof(bool)* MAX_BONDS);
	sbs = (SimBond *)malloc(sizeof(SimBond)* MAX_BONDS*MAX_BONDS);
	//the_mat = (int *)malloc(sizeof(int *)*MAX_BONDS*MAX_BONDS);
	the_mat =(int **)malloc(sizeof(int *)*MAX_BONDS);
	for(int i=0; i<MAX_BONDS; i++) the_mat[i]=(int *)malloc(sizeof(int)*MAX_BONDS);




	blev_mat = calcBondLevelsMatrix(parms);
	mi1=(MolInfo *)malloc(sizeof(MolInfo));
	mi2=(MolInfo *)malloc(sizeof(MolInfo));


	if((fil1 = fopen( fname1, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname1 );
		exit(1);
	}
	eof=readMolInfo(fil1, mi1);
	natoms1 = mi1->natoms;
	fclose(fil1);

	//initialize the init weights - if no weigths are available then initialize them all at 1.0
	if(iwtsFname[0]) {
		init_wts = readAtomWeights(iwtsFname, natoms1); 
	} else {
		init_wts=(float *)malloc(sizeof(float) * natoms1);
		for(int row=0; row<natoms1; row++) init_wts[row]=1.0f;
	}
	ss1 = getDiagonalElements_SS(mi1, parms, blev_mat);

	ss2	= (float *)malloc(sizeof(float) * num_mols);
	natoms2 = (int *)malloc(sizeof(int) * num_mols); //number of atoms for each molecule
	matrices = (float ***)malloc(sizeof(float **) * num_mols);


	if((fil2 = fopen( fname2, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname2 );
		exit(1);
	}

	//here it starts 
	eof=readMolInfo(fil1, mi2);
	m=0;
	while(!eof) {
		goforit=true;
		if(filt_db) goforit=cid_in_flist(mi2->cid, filt_db);
		if(goforit==true) {
			//save the number of atoms
			natoms2[m] = mi2->natoms;
			//printf("%d  --> %d\n", mi2->cid, mi2->natoms); 
			//get also the maximum possible number of atoms in any molecule
			if(mi2->natoms > max_natoms2) max_natoms2=mi2->natoms;
			//puts("4a");
			//compute the self similarities
			//printf("%d  --> %d\n", mi2->cid, mi2->natoms);
			ss2[m]=calcSelfSimilarity(mi2, parms, blev_mat)/10000.0f;
			//puts("4b");
			//get the aba_match_scores for the 2 mols
			aba_match_scores = getAbaBondsCompMatrix(mi1, mi2, parms); 
			//puts("4c");
			//compute the atom matrix <- this is the important stuff to be saved
			atom_matrix=calcAtomMatchingMatrix(mi1, mi2, aba_match_scores, blev_mat, parms->PEN);
			matrices[m]=array_to_matrix_float(atom_matrix, mi1->natoms, mi2->natoms);
			free(aba_match_scores); //clean up
			free(atom_matrix);
			m++;
		}
		freeMolInfo(mi2);
		eof=readMolInfo(fil2, mi2);
	}
	fclose(fil2);

	//now the cruel part
	//first allocate space for the weights
	wts=(float *)malloc(sizeof(float)*natoms1);
	//alocate space also for the scores
	scores=(float *)malloc(sizeof(float)*m);
	//allocate the space to save the work matrix using as many cols as the largest possible number of atoms in themolecule list
	int** mat=(int **)malloc(sizeof(int*)*natoms1);
	for(int row=0; row<natoms1; row++) mat[row]=(int *)malloc(sizeof(int) * max_natoms2);
	srand( (unsigned)time( NULL ) );
	for(int i=0;i<nruns; i++){
		//generate the weights
		sumw=0;
		for(int a=0;a<natoms1; a++) {
			rnd=((float)rand()) / (RAND_MAX + 1); //values between 0 and 1.0
			if(i>0) {
				wts[a]=fabsf(init_wts[a] + 2.0f*rnd*spread - spread); //no negative weights ever!
				//if(wts[a]>2.0f) wts[a]=2.0f;
			} else wts[a]=init_wts[a];  //the null option
			//if(wFile) fprintf(wFile, "%7.4f ", wts[a]);
			sumw+=wts[a];
		}
		fac=((float)natoms1)/sumw;
		for(int a=0;a<natoms1; a++) {
			if(i>0) wts[a]*=fac;
			if(wFile) fprintf(wFile, "%7.4f ", wts[a]);

		}

		if(wFile) fprintf(wFile,"\n");

		WeightsEvaluator(matrices, mat, wts, m, natoms1, natoms2,  ss1, ss2, scores, J);
		for(int mol=0; mol<m; mol++) {
			printf("%7.4f ",scores[mol]);
		}
		printf("\n");
	}

	free(wts);
	free(scores);
	free_array_to_matrix(mat, natoms1); 
	freeMolInfo(mi1); 
	free(mi1);
	free(mi2);

	free(blev_mat);
	return 1;
	
}
