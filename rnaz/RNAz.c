/*********************************************************************                
 *                                                                   *
 *                              RNAz.c                               *
 *                                                                   *
 *	Assess alignments for exceptionally stable and/or conserved  *
 *	secondary structures using RNAalifold/RNAfold and SVMs.      *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                  *
 *                                                                   *
 *	   $Id: RNAz.c,v 1.5 2004-10-21 08:02:04 wash Exp $          *
 *                                                                   *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
#include "alifold.h"
#include "zscore.h"
#include "svm.h"
#include "svm_helper.h"

#define PRIVATE static
#define MAX_NUM_NAMES 500
#define MIN2(A, B) ((A) < (B) ? (A) : (B))

PRIVATE void usage(void);
PRIVATE void help(void);
PRIVATE void version(void);

PRIVATE int read_clustal(FILE *clust,
						 char *alignedSeqs[],char *names[]);

PRIVATE char *consensus(const char *AS[]);

PRIVATE double meanPairID(const char *AS[]);

PRIVATE void revAln(char *AS[]);

PRIVATE void sliceAln(const char *sourceAln[], char *destAln[],
					  int from, int to);

PRIVATE void freeAln(char *AS[]);

PRIVATE void classify(double* prob, double* decValue,
					  struct svm_model* decision_model,
					  double id,int n_seq, double z,double sci);




/********************************************************************
 *                                                                  *
 * main -- main program                                             *
 *                                                                  *
 ********************************************************************/

int main(int argc, char *argv[])
{

  char *modelDir=NULL;              /* Directory with model files */
  struct svm_model* decision_model; /* SVM classification model */

  /* Command line options */
  int reverse=0;     /* Scan reverse complement */
  int showVersion=0; /* Shows version and exits */
  int showHelp=0;    /* Show short help and exits */
  int from=-1;       /* Scan slice from-to  */
  int to=-1;

  FILE *clust_file=stdin; /* Input file */

  /* Arrays storing sequences/names of the alignments */
  char *AS[MAX_NUM_NAMES];     
  char *names[MAX_NUM_NAMES];
  char *window[MAX_NUM_NAMES]; 

  int n_seq;  /* number of input sequences */
  int length; /* length of alignment/window */

  char *structure=NULL;
  char *singleStruc, *output,*woGapsSeq;
  char strand[8];
  char *string=NULL;
  double singleMFE,sumMFE,singleZ,sumZ,z,sci,id,decValue,prob;
  double min_en, real_en;
  int i,j,r,countSeq;

  /* Global RNA package variables */
  do_backtrack = 1; 
  dangles=2;

  
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') {
      switch ( argv[i][1] )	{
	  case 'r': reverse=1;break;
	  case 'f':
		r=sscanf(argv[++i], "%d", &from);
		if (r!=1) usage();
		break;
	  case 't':
		r=sscanf(argv[++i], "%d", &to);
		if (r!=1) usage();
		break;
	  case 'v':showVersion=1;break;
	  case 'h':showHelp=1;break;
	  default: usage();
	  }
    }
    else { /* doesn't start with '-' should be filename */ 
      if (i!=argc-1) usage();
      clust_file = fopen(argv[i], "r");
      if (clust_file == NULL) {
		fprintf(stderr, "ERROR: Can't open %s\n", argv[i]);
		usage();
      }
    }
  }

  if (showVersion==1) version();
  if (showHelp==1) help();

  
  modelDir=getenv("RNAZDIR");

  if (modelDir==NULL){
	nrerror("ERROR: Could not find models. You have to set the RNAZDIR"
			" enviroment variable pointing to the model files!\n");
  }
 
  //n_seq = read_clustal(clust_file, AS, names);

  //if (clust_file != stdin) fclose(clust_file);
  //if (n_seq==0) nrerror("ERROR: There were problems"
  //					" reading the input file.\n");


  countSeq=0;

  regression_svm_init(modelDir);

  decision_model=get_decision_model(modelDir);

  if (decision_model==NULL){
	nrerror("ERROR: Could not find decision-model. You have to set "
			"the RNAZDIR enviroment variable pointing to the model files!\n");
  }


  
  while ((n_seq=read_clustal(clust_file, AS, names))!=0){
	countSeq++;
	length = (int) strlen(AS[0]);
  
	/* if a slice is specified by the user */
  
	if ((from!=-1 || to!=-1) && (countSeq==1)){
	  if (((from!=-1) && (to==-1)) || ((from==-1) && (to!=-1))){
		nrerror("ERROR: You have to set both -f and -t parameters"
				" to score a slice of the alignment.\n");
	  }
	  if ((from>=to)||(from<=0)||(to>length)){
		nrerror("ERROR: -f and/or -t parameters out of range"
				" (no reasonable slice specified)\n");
	  }
	  sliceAln((const char **)AS, (char **)window, from, to);
	  length=to-from+1;
	} else { /* take complete alignment */
	  /* window=AS does not work..., deep copy seems not necessary here*/
	  from=1;
	  to=length;
	  sliceAln((const char **)AS, (char **)window, 1, length);
	}
  
	if (reverse==1){
	  revAln((char **)window);
	  strcpy(strand,"reverse");
	} else {
	  strcpy(strand,"forward");
	}
  
	structure = (char *) space((unsigned) length+1);
	min_en = alifold(window, structure);

	sumZ=0;
	sumMFE=0;

	output=(char *)space(sizeof(char)*(length+16)*(n_seq+1)*3);
  
	for (i=0;i<n_seq;i++){
	  singleStruc = space(strlen(window[i])+1);
	  woGapsSeq = space(strlen(window[i])+1);
	  j=0;
	  while (window[i][j]){
		window[i][j]=toupper(window[i][j]);
		if (window[i][j]=='T') window[i][j]='U';
		if (window[i][j]!='-'){
		  woGapsSeq[strlen(woGapsSeq)]=window[i][j];
		  woGapsSeq[strlen(woGapsSeq)]='\0';
		}
		++j;
	  }
	  
	  singleMFE = fold(woGapsSeq, singleStruc);
	  singleZ=mfe_zscore(woGapsSeq,singleMFE);
	  
	  sumZ+=singleZ;
	  sumMFE+=singleMFE;
	  sprintf(output+strlen(output),">%s\n%s\n%s ( %6.2f)\n",
			  names[i],woGapsSeq,singleStruc,singleMFE);
	  free(woGapsSeq);
	  free(singleStruc);
	}

	{
	  int i; double s=0;
	  extern int eos_debug;
	  eos_debug=-1; /* shut off warnings about nonstandard pairs */
	  for (i=0; window[i]!=NULL; i++) 
		s += energy_of_struct(window[i], structure);
	  real_en = s/i;
	}

	string = consensus((const char **) window);
	sprintf(output+strlen(output),
			">consensus\n%s\n%s (%6.2f = %6.2f + %6.2f) \n",
			string, structure, min_en, real_en, min_en-real_en );

	id=meanPairID((const char **)window);
	z=sumZ/n_seq;
	sci=min_en/(sumMFE/n_seq);
	
	decValue=999;
	classify(&prob,&decValue,decision_model,id,n_seq,z,sci);
	printf("\n###########################  RNAz 0.1.1  #############################\n\n");
	printf(" Sequences: %u\n", n_seq);
	printf(" Slice: %u to %u\n",from,to);
	printf(" Columns: %u\n",length);
	printf(" Strand: %s\n",strand);
	printf(" Mean pairwise identity: %6.2f\n", id);
	printf(" Mean single sequence MFE: %6.2f\n", sumMFE/n_seq);
	printf(" Consensus MFE: %6.2f\n",min_en);
	printf(" Energy contribution: %6.2f\n",real_en);
	printf(" Covariance contribution: %6.2f\n",min_en-real_en);
	printf(" Mean z-score: %6.2f\n",z);
	printf(" Structure conservation index: %6.2f\n",sci);
	printf(" SVM decision value: %6.2f\n",decValue);
	printf(" SVM RNA-class probability: %6f\n",prob);
	if (prob>0.5){
	  printf(" Prediction: RNA\n");
	}
	else {
	  printf(" Prediction: no RNA\n");
	}
	
	printf("\n######################################################################\n\n");
	
	printf("%s//\n\n",output);
	
	fflush(stdout);
	
	free(output);
	freeAln((char **)AS);
	freeAln((char **)names);
	freeAln((char **)window);
	

  }

  svm_destroy_model(decision_model);
  regression_svm_free();
  
  return 0;
}

/********************************************************************
 *                                                                  *
 * read_clustal -- read CLUSTAL W formatted file                    *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * clust ... filehandle pointing to the file to be read             *
 * alignedSeqs ... array of strings where read sequences are stored *
 * names ... array of sequence names                                *
 *                                                                  *
 * Returns number of sequences read                                 *
 *                                                                  *
 ********************************************************************/


PRIVATE int read_clustal(FILE *clust, char *alignedSeqs[], char *names[]) {

  char *line, name[100]={'\0'}, *seq;
  int  n, nn=0, num_seq = 0;
   
  if ((line=get_line(clust)) == NULL) {
	//fprintf(stderr, "ERROR: Empty CLUSTAL file\n");
	return 0;
  }

  
  while (((n=strlen(line))<4) || isspace((int)line[0])){
	free(line); line = get_line(clust);
  }

 
  if (strncmp(line,"CLUSTAL", 7) !=0) {
	fprintf(stderr, "ERROR: No CLUSTAL file\n");
	free(line); return 0;
  }
  free(line);
  line = get_line(clust);

  while ((line!=NULL) && (strcmp(line,"//")!=0)) {
	if (((n=strlen(line))<4) || isspace((int)line[0])) {
	  /* skip non-sequence line */
	  free(line); line = get_line(clust);
	  nn=0; /* reset seqence number */
	  continue;
	} 
     
	seq = (char *) space( (n+1)*sizeof(char) );
	sscanf(line,"%99s %s", name, seq);
	if (nn == num_seq) { /* first time */
	  names[nn] = strdup(name);
	  alignedSeqs[nn] = strdup(seq);
	}
	else {
	  if (strcmp(name, names[nn])!=0) {
		/* name doesn't match */
		fprintf(stderr,
				"ERROR: Inconsitent sequence names in CLUSTAL file)\n");
		free(line); free(seq);
		return 0;
	  }
	  alignedSeqs[nn] = (char *)
		xrealloc(alignedSeqs[nn], strlen(seq)+strlen(alignedSeqs[nn])+1);
	  strcat(alignedSeqs[nn], seq);
	}
	nn++;
	if (nn>num_seq) num_seq = nn;
	free(seq);
	free(line);
	if (num_seq>=MAX_NUM_NAMES) {
	  fprintf(stderr, "ERROR: Too many sequences in CLUSTAL file\n");
	  return 0;
	}
	line = get_line(clust);
  }

  alignedSeqs[num_seq] = NULL;
  names[num_seq] = NULL;
  if (num_seq == 0) {
	fprintf(stderr, "ERROR: No sequences found in CLUSTAL file\n");
	return 0;
  }
  n = strlen(alignedSeqs[0]); 
  for (nn=1; nn<num_seq; nn++) {
	if (strlen(alignedSeqs[nn])!=n) {
	  fprintf(stderr, "ERROR: Sequences are of unequal length.\n");
	  return 0;
	}
  }
  return num_seq;
}
 

/********************************************************************
 *                                                                  *
 * consensus -- Calculates consensus of alignment                   *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * AS ... array with sequences                                      *
 *                                                                  *
 * Returns string with consensus sequence                           *
 *                                                                  *
 ********************************************************************/

PRIVATE char *consensus(const char *AS[]) {
  char *string;
  int i,n;
  n = strlen(AS[0]);
  string = (char *) space((n+1)*sizeof(char));
  for (i=0; i<n; i++) {
    int s,c,fm, freq[8] = {0,0,0,0,0,0,0,0};
    for (s=0; AS[s]!=NULL; s++) 
      freq[encode_char(AS[s][i])]++;
    for (s=c=fm=0; s<8; s++) /* find the most frequent char */
      if (freq[s]>fm) {c=s, fm=freq[c];}
    if (s>4) s++; /* skip T */
    string[i]=Law_and_Order[c];
  }
  return string;
}


/********************************************************************
 *                                                                  *
 * meanPairID -- Calculates mean pairwise identity of alignment     *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * AS ... array with sequences                                      *
 *                                                                  *
 * Returns mean pair ID in percent                                  *
 *                                                                  *
 ********************************************************************/


PRIVATE double meanPairID(const char *AS[]) {

  int i,j,k,matches,pairs,length;

  matches=0;
  pairs=0;

  length=strlen(AS[0]);
  
  for (i=0;AS[i]!=NULL;i++){
	for (j=i+1;AS[j]!=NULL;j++){
	  for (k=0;k<length;k++){
		if ((AS[i][k]!='-') || (AS[j][k]!='-')){
		  if (AS[i][k]==AS[j][k]){
			matches++;
		  }
		  pairs++;
		}
	  }
	}
  }

  return (double)(matches)/pairs*100;
  
}

/********************************************************************
 *                                                                  *
 * revAln -- Reverse complements sequences in an alignment          *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * AS ... array with sequences which is rev-complemented in place   *
 *                                                                  *
 ********************************************************************/

void revAln(char *AS[]) {

  int i,j,length;
  char *tmp;
  char letter;
  length=strlen(AS[0]);
  
  for (i=0;AS[i]!=NULL;i++){
	tmp = (char *) space((unsigned) length+1);
	for (j=length-1;j>=0;j--){
	  letter=AS[i][j];
	  switch(letter){
		case 'T': letter='A'; break;
		case 'U': letter='A'; break;
		case 'C': letter='G'; break;
		case 'G': letter='C'; break;
		case 'A': letter='U'; break;
	  }
	  tmp[length-j-1]=letter;
	}
	tmp[length]='\0';
	strcpy(AS[i],tmp);
	free(tmp);
	tmp=NULL;
  }
}

/********************************************************************
 *                                                                  *
 * sliceAln -- Gets a slice of an alignment                         *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * sourceAln ... array with sequences of source alignment           *
 * destAln   ... pointer to array where slice is stored             *
 * from, to  ... specifies slice, first column is column 1          *
 *                                                                  *
 ********************************************************************/

void sliceAln(const char *sourceAln[], char *destAln[],
			  int from, int to){

  int i;
  char *slice;

  for (i=0;sourceAln[i]!=NULL;i++){
	slice=(char *) space((unsigned) (to-from+2));
	strncpy(slice,sourceAln[i]+from-1,(to-from+1));
	destAln[i]=slice;
  }
  destAln[i]=NULL;
}

/********************************************************************
 *                                                                  *
 * freeAln -- Frees memory of alignment array                       *
 *                                                                  *
 ********************************************************************/

PRIVATE void freeAln(char *AS[]){

  int i;
  for (i=0;AS[i]!=NULL;i++){
	free(AS[i]);
  }
  AS=NULL;
}



/********************************************************************
 *                                                                  *
 * classify -- SVM classification depending on various variables    *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * prob ... pointer where class probability is stored               *
 * decValue ... pointer where decidion value is stored              *
 * svm_model ... pointer to SVMLIB model used for the classificaton *
 * id ... mean pairwise identity of alignment                       *
 * n_seq ... number of sequences in the alignment                   *
 * z ... mean z-score of single sequences in the alignment          *
 * sci ... structure conservation index of alignment                *
 *                                                                  *
 ********************************************************************/


PRIVATE void classify(double* prob, double* decValue,
					  struct svm_model* decision_model,
					  double id,int n_seq, double z,double sci){

  struct svm_node node[5];

  double value=0;
  
  node[0].index = 1; node[0].value = z;
  node[1].index = 2; node[1].value = sci;
  node[2].index = 3; node[2].value = id;
  node[3].index = 4; node[3].value = n_seq;
  node[4].index =-1;

  scale_decision_node((struct svm_node*)&node);

  svm_predict_values(decision_model,node,&value);
  *decValue=value;
  svm_predict_probability(decision_model,node,&value);
  *prob=value;
}


/********************************************************************
 *                                                                  *
 * usage, help, version - shows information and exits               *
 *                                                                  *
 ********************************************************************/

PRIVATE void usage(void)
{
  nrerror("RNAz [-v] [-h] [-r] [-f from] [-t to] file\n\n"
		  "Run RNAz -h for details.\n");
}

PRIVATE void help(void){

  nrerror("RNAz [-v] [-h] [-r] [-f from] [-t to] file\n\n"
		  "file    Input alignment in Clustal W format\n"
		  "-f -t   Score subregion from-to of alignment\n"
		  "-r      Scan reverse complement of input alignment\n"
		  "-v      Show version information\n"
		  "-h      Show this help message\n"
		  );

}

PRIVATE void version(void){
  printf("RNAz 0.1.1, September 2004\n");
  exit(EXIT_SUCCESS);
}

