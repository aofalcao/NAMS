
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "nams.h"
#include "awts.h"



#define MOL_AGAINST_MOL 0
#define MOL_AGAINST_LIST 1
#define LIST_AGAINST_ITSELF 2
#define LIST_AGAINST_LIST 3
#define MOL_AGAINST_LIST_W 4

bool USE_MW=true;
float THRES_JACC=0.0f;



bool charIn(char c, char *targ) 
{
	unsigned int i=0;
	while(i<strlen(targ)) if(c==targ[i++]) return true;
	return false;
}

int displayHelp() 
{

puts("NAMS - Non-contiguous atom matching structural similarity function");
printf("nams version - %s\n", NAMS_VERSION);
char s1[]="This program will compare one or several molecules against others producing\n\
a similaity score for each pair of molecules\n\
Basic use:\n\
nams -mode [mode] -in [mol_filename] -db [mod db filename]\n\
NOTE: -mode, and -db are required options\n\
Details\n\
    -in [filename] filename is the name of a file in the NAMS format with one\n\
       or more molecules   \n\
    -db [filename] filename is the name of a file in the NAMS format with one\n\
        or more molecules   \n\
    -mode: Operation mode - must follow a value: 1,2,3 or 4\n\
        1 - Compares one molecule (given with the -in option) with a full \n\
            database (-db)\n\
        2 - Compares a database with itself. All molecules of db are compared\n\
            to all the others in the database\n\
        3 - Compares a list of molecules (-in) with all the molecules in (-db)";
char s2[]="    -opts [AJLMS] Output options. Each of these 5 letters or combinations may \n\
            be used. If options A or M are absent results are provided in a table\n\
            [default -JLS - if -opts present all is reset]\n\
        M - Presents the atom matching scores for all atom pairs for both \n\
            molecules\n\
        A - Gives the optimum matching between atoms of both molecules. \n\
            Each atom pair is identified by their atom ID and the \n\
            corresponding matching score\n\
        L - Self similarity scores of both molecules [default]\n\
        S - Similarity score between both molecules [default]\n\
        J - Jaccard score between both molecules [default]\n\
Filter parameters:\n\
    -tj [value] - presents only the results above a given Jaccard threshold\n\
    -fdb [filename] -defines a filename that includes the molids to search within\n\
	     the database, as specified by -db \n\
    -fin [filename] -defines a filename that includes the molids to search within\n\
	     the NAMS input file, as specified by -in\n\
    -fmxmw [value] - the maximum molecular %weight difference to search within the \n\
             database [default=9999]\n\
    -fmnmw [value] - the minimum molecular %weight difference to search within the \n\
             database [default=100]\n\
NAMS Heuristic parameters:\n\
    -alpha [value] - alpha parameter of the model [default=0.9]\n\
    -nrings [value] - score that acounts for how many rings an atom belongs\n\
           [default=0.80]\n\
    -stiso [value] - parameter for accountin chiral stereo isomerism\n\
           [default=0.95]\n\
    -ctiso [value] - accounts for double bond cis-trans isomerisms\n\
           [default=0.95]\n\
    -bring [value] - whether a bond is in a ring or not [default=0.90]\n\
    -barom [value] - whether a bond is in an aromatic ring or not\n\
           [default = 0.90]\n\
    -border [value] - bond covalent order [default=0.90]\n\
    -pen [value] - penalty factor to account for unmatched bonds \n\
           [default= 1.0]\n\
    -ADM [code] - atom distance matrix code. Admissible values are 0,1,2, 3 and 4\n\
           [default=3]"; 

    char s3[]="(c) Andre Falcao 2013-2016 \n\
    This Program is provided free of charge and WITHOUT ANY WARRANTY\n\
  if used please cite:\n\
     AL Teixeira, AO Falcao. 2013. Non contiguous atom matching structural similarity\n\
     function DOI: 10.1021/ci400324u. Journal of Chemical Information and Modeling";
char o4[]= "        4 - Compares a molecule against a database, by randomly generating a set of\n\
            atom weights to the base molecule. Outputs only the Jaccard results for\n\
            each molecule in a column with as many rows as the number of runs";

char sw[]="Atom Weights: (modes 0,1 and 4)\n\
    -iwfname [name] - name of the input file with sigle weights for each atom of the\n\
	       input molecule [default = all weights are 1.0]\n\
    -owfname [name] - name of output file with the generated weights\n\
    -spread [value] - spread around the default values from which new values will be\n\
	       randommly generated [default=0.2]\n\
    -nruns [value] - number of runs for atomic weights generation [default=100]";


	puts(s1);
	puts(o4);	//to remove
	puts(s2);
	puts(sw);	//to remove
	puts(s3);

	return 1;
} 

int *ReadFilteredMolecules(char *fname) {
	//reads the molids to be considered in a search - file must have one molid per row
	//returns 
	FILE *fil;
	char line[MAX_LIN_SIZ];
	int count=0;
	int *flist;
	if((fil = fopen( fname, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname );
		exit(1);
	}
	//I will read it twice, to get all the molids first and then create an array and fill it
	while(!feof(fil)) {
		fgets(line, MAX_LIN_SIZ,fil);
		count++;
	}
	fclose(fil);
	if((fil = fopen( fname, "rt" )) == NULL ) {
		printf( "The file '%s' was not opened - It's poorly formated or missing!\n", fname );
		exit(1);
	}
	flist=(int *)malloc((count+1)*sizeof(int));
	//the first item is the COUNT of molecules read
	flist[0]=count;
	count=1;
	while(!feof(fil)) {
		fgets(line, MAX_LIN_SIZ,fil);
		flist[count++]=atoi(line);
	}
	fclose(fil);
	//printf("---->%d\n", flist[0]);
	return flist;
}



int main(int argc, char *argv[]) 
{
    clock_t start, end;
	Params  parms;
	int		p=0, error_flag=0, ADM=3;
	char	fname1[MAX_LIN_SIZ], fname2[MAX_LIN_SIZ];
	char	fname_filtin[MAX_LIN_SIZ], fname_filtdb[MAX_LIN_SIZ];
	int		*filt_in=NULL, *filt_db=NULL;
	int		min_mol_wt=100, max_mol_wt=9999;
	float	fmin_mol_wt, fmax_mol_wt;
	char	opts[8];
	int		mode=MOL_AGAINST_MOL;
	bool	bin=false, bdb=false, bopts=false, silence;
	float	BS_ALPHA = 0.9f;
    float	ANRINGS_FAC =0.8f;
    float	ACHIR_FAC = 0.95f;    
    float	DBSTEREO_FAC = 0.95f;
    float	BRING_FAC = 0.9f;     
    float	BAROM_FAC = 0.9f;     
    float	BORDER_FAC = 0.9f;    
    float	PEN = 1.0f;
	char	iwtsFname[MAX_LIN_SIZ];
	char	owtsFname[MAX_LIN_SIZ];
	float	spread=0.1f;
	owtsFname[0]='\0';
	iwtsFname[0]='\0';

	int nruns=100;
	int palpha=0, pnrings=0, pstiso=0, pctiso=0, pbring=0, pbarom=0, 
		pborder=0, ppen=0, pADM=0, pin=0, pdb=0, popts=0, pmode=0, amode,
		pfmxmw=0, pfmnmw=0, pfdb=0, pfin=0, ptj=0, pnruns=0, powfname=0,
		piwfname=0, pspread=0;
	bool S=true, L=true, J=true, M=false, A=false;
	//puts("hi there!");
	for(int p=1; p < argc; p++) {
		if(!strcmp(argv[p],"-mode")) pmode=p+1;
		if(!strcmp(argv[p],"-alpha")) palpha=p+1;

		if(!strcmp(argv[p],"-fdb")) pfdb=p+1;
		if(!strcmp(argv[p],"-fin")) pfin=p+1;
		if(!strcmp(argv[p],"-fmxmw")) pfmxmw=p+1;
		if(!strcmp(argv[p],"-fmnmw")) pfmnmw=p+1;
		 
		if(!strcmp(argv[p], "-nrings")) pnrings=p+1;
		if(!strcmp(argv[p], "-stiso")) pstiso=p+1;
		if(!strcmp(argv[p], "-ctiso")) pctiso=p+1;
		if(!strcmp(argv[p], "-bring")) pbring=p+1;
		if(!strcmp(argv[p], "-barom")) pbarom=p+1;
		if(!strcmp(argv[p], "-border")) pborder=p+1;
		if(!strcmp(argv[p], "-pen")) ppen=p+1;
		if(!strcmp(argv[p], "-ADM")) pADM=p+1;
		if(!strcmp(argv[p], "-in")) pin=p+1;
		if(!strcmp(argv[p], "-db")) pdb=p+1;
		if(!strcmp(argv[p], "-opts")) popts=p+1;
		if(!strcmp(argv[p], "-tj")) ptj=p+1; 
		if(!strcmp(argv[p], "-nruns")) pnruns=p+1; 
		if(!strcmp(argv[p], "-owfname")) powfname=p+1; 
		if(!strcmp(argv[p], "-iwfname")) piwfname=p+1; 
		if(!strcmp(argv[p], "-spread")) pspread=p+1; 
		if(!strcmp(argv[p], "-OLD")) USE_MW=false; //old format does not include mol weight 
		if(!strcmp(argv[p], "-silent")) silence=true;
		if(!strcmp(argv[p], "-help")) {
			displayHelp();
			return 0;
		}
		
	}
	  
	if(pmode>0) {
		amode=atoi(argv[pmode]);
		switch(amode) {
			case 1: 
				mode = MOL_AGAINST_LIST;
				break;
			case 2:
				mode = LIST_AGAINST_ITSELF;
				break;
			case 3:
				mode = LIST_AGAINST_LIST;
				break;
			case 4:   //to remove
				mode = MOL_AGAINST_LIST_W;
				break;
		}
	}
	if(pin > 0) {
		bin=true;
		strcpy(fname1, argv[pin]);
	}
	if(pdb > 0) {
		bdb=true;
		strcpy(fname2, argv[pdb]);
	}
	if(popts > 0) {
		bopts=true;
		S=false;
		L=false;
		J=false;
		M=false;
		A=false;
		strcpy(opts, argv[popts]);
		//a) their similarity scores (S)
		//b) their self similarity scores (L)
		//c) their Jaccard coefficient (J)
		//d) the atom matching matrix (M)
		//e) the alignment produced (A)
		if(charIn('S', opts)) S=true;
		if(charIn('L', opts)) L=true;
		if(charIn('J', opts)) J=true;
		if(charIn('M', opts)) M=true;
		if(charIn('A', opts)) A=true;
	}
	if(mode==LIST_AGAINST_ITSELF) {
		if( !bdb) {
			puts("Missing  -db - This option is required!");
			puts("See nams -help  -  for instructions of use");
			exit(1);
		}
	} else {
		if(!bin || !bdb) {
			puts("Missing  -in or -db - These options are required!");
			puts("See nams -help  -  for instructions of use");
			exit(1);
		}
	}
	if(ptj > 0) THRES_JACC=(float)atof(argv[ptj]);
	if(palpha > 0) BS_ALPHA=(float)atof(argv[palpha]);
	if(pnrings > 0) ANRINGS_FAC=(float)atof(argv[pnrings]);
	if(pstiso > 0) ACHIR_FAC=(float)atof(argv[pstiso]);
	if(pctiso > 0) DBSTEREO_FAC=(float)atof(argv[pctiso]);
	if(pbring > 0) BRING_FAC=(float)atof(argv[pbring]);
	if(pbarom > 0) BAROM_FAC=(float)atof(argv[pbarom]);
	if(pborder > 0) BORDER_FAC=(float)atof(argv[pborder]);
	if(ppen > 0) PEN=(float)atof(argv[ppen]);
	if(pADM > 0) ADM=atoi(argv[pADM]);

	if(pfin > 0) {
		strcpy(fname_filtin, argv[pfin]);
		filt_in=ReadFilteredMolecules(fname_filtin);
	}
	if(pfdb > 0) {
		strcpy(fname_filtdb, argv[pfdb]);
		filt_db=ReadFilteredMolecules(fname_filtdb);
	}
	//printf("%d %d %s %s\n", min_mol_wt, max_mol_wt, fname_filtin,fname_filtdb);

	if(pfmxmw > 0) max_mol_wt=atoi(argv[pfmxmw]);
	if(pfmnmw > 0) min_mol_wt=atoi(argv[pfmnmw]);
	if(pnruns > 0) nruns=atoi(argv[pnruns]);
	if(powfname > 0) strcpy(owtsFname, argv[powfname]);
	if(piwfname > 0) strcpy(iwtsFname, argv[piwfname]);
	if(pspread > 0) spread=(float)atof(argv[pspread]);
	//printf("------> %f\n", spread);

	setParams(&parms, ADM, BS_ALPHA, ANRINGS_FAC, ACHIR_FAC, DBSTEREO_FAC, BRING_FAC, BAROM_FAC, BORDER_FAC, PEN);

    start = clock();
	fmin_mol_wt=(100.0f-(float)min_mol_wt)/100.0f;
	fmax_mol_wt=(100.0f+(float)max_mol_wt)/100.0f;
	//printf("%f %f\n", fmin_mol_wt, fmax_mol_wt);
	if(mode == MOL_AGAINST_MOL) Simple2Molecules(fname1, fname2, &parms, S,L,J,M,A, iwtsFname);
	if(mode == MOL_AGAINST_LIST) MoleculeAgainstList(fname1, fname2, &parms, S,L,J,M,A, fmin_mol_wt, fmax_mol_wt, filt_db, iwtsFname, silence);
	if(mode == LIST_AGAINST_ITSELF) ListAgainstItself(fname2, &parms, S,L,J,M,A, fmin_mol_wt, fmax_mol_wt, filt_in, filt_db, silence); 
	if(mode == LIST_AGAINST_LIST) ListAgainstList(fname1, fname2, &parms, S,L,J,M,A, fmin_mol_wt, fmax_mol_wt, filt_in, filt_db,silence); 
	if(mode == MOL_AGAINST_LIST_W) GenerateAtomWeights(fname1, fname2, &parms, J, filt_db, nruns, owtsFname, iwtsFname, spread); 

	
	end = clock();
    if(!silence) fprintf(stderr, "Time required for execution: %f seconds\n",  (double)(end-start)/CLOCKS_PER_SEC);
	return 0;	
}

