/* Last changed Time-stamp: <04/09/17 17:45:30 wash> */
/*                
	Assess alignments for exceptionally stable and/or conserved
	secondary structures using alifold-routines and SVMs.

	               c Stefan Washietl, Ivo L Hofacker
              		  
*/

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

static const char rcsid[] = "$Id: RNAz.c,v 1.1 2004-09-18 13:25:52 wash Exp $";

#define PRIVATE static

PRIVATE void usage(void);
PRIVATE void help(void);
PRIVATE int read_clustal(FILE *clust, char /*@out@*/ *AlignedSeqs[],
			 char /*@out@*/ *names[]);
PRIVATE char *consensus(const char *AS[]);
PRIVATE double meanPairID(const char *AS[],int n_seq);
PRIVATE void revAln(char *AS[],int n_seq);
PRIVATE void sliceAln(const char *sourceAln[], char *destAln[], int n_seq, int from, int to);
PRIVATE void classify(double* prob, double* decValue,struct svm_model* decision_model,double id,int n_seq, double z,double sci);

#define MAX_NUM_NAMES    500

int main(int argc, char *argv[])
{
  char *string=NULL;
  char *structure=NULL;
  char *singleStruc, *output,*woGapsSeq;
  char *modelDir=NULL;
  int j;
  double singleMFE,sumMFE,singleZ,sumZ,z,sci,id,decValue,prob;
  struct svm_model* decision_model;
  int   n_seq, i, length, from, to,r;
  double min_en, real_en;
  int   istty;
  char     *AS[MAX_NUM_NAMES];          /* aligned sequences */
  char     *names[MAX_NUM_NAMES];       /* sequence names */
  char     *window[MAX_NUM_NAMES];
  FILE     *clust_file = stdin;
  char strand[8];
  int reverse=0, showVersion=0, showHelp=0;
  
  do_backtrack = 1; 
  dangles=2;
  from=-1; to=-1;
  
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
		fprintf(stderr, "can't open %s\n", argv[i]);
		usage();
      }
    }
  }

  if (showVersion==1){
	printf("RNAz 0.1.1, September 2004\n");
	exit(0);
  }

  if (showHelp==1) help();

  
  modelDir=getenv("RNAZDIR");

  if (modelDir==NULL){
	nrerror("RNAz 0.1 Error:\nCould not find models. You have to set the RNAZDIR enviroment variable pointing to the model files!\n");
  }
  
  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
  
  n_seq = read_clustal(clust_file, AS, names);

  if (clust_file != stdin) fclose(clust_file);
  if (n_seq==0)
    nrerror("RNAz 0.1 Error:\nNo sequences found.\n");

  
  length = (int) strlen(AS[0]);

  /* if a slice is specified by the user */
  if (from!=-1 || to!=-1){
  	if ((from!=-1) && (to==-1) || (from==-1) && (to!=-1)){
	  nrerror("RNAz 0.1 Error:\nYou have to set both -f and -t parameters to score a slice of the alignment.\n");
	}
	if ((from>=to)||(from<=0)||(to>length)){
	  nrerror("RNAz 0.1 Error:\n -f and/or -t parameters out of range (no reasonable slice specified)\n");
	}
	sliceAln((const char **)AS, (char **)window, n_seq, from, to);
	length=to-from+1;
  } else { /* take complete alignment */
	/* window=AS does not work..., deep copy seems not necessary here*/
	from=1;
	to=length;
	sliceAln((const char **)AS, (char **)window, n_seq, from, length);
  }
  
  if (reverse==1){
	revAln((char **)window,n_seq);
	strcpy(strand,"reverse");
  } else {
	strcpy(strand,"forward");
  }
  
  id=meanPairID((const char **)window,n_seq);

  structure = (char *) space((unsigned) length+1);
  min_en = alifold(window, structure);

  sumZ=0;
  sumMFE=0;

  output=(char *)space(sizeof(char)*(length+16)*(n_seq+1)*3);
  
  regression_svm_init(modelDir);

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
	sprintf(output+strlen(output),">%s\n%s\n%s ( %6.2f)\n",names[i],woGapsSeq,singleStruc,singleMFE);
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
  sprintf(output+strlen(output),">consensus\n%s\n%s (%6.2f = %6.2f + %6.2f) \n",string, structure, min_en, real_en, min_en-real_en );

  z=sumZ/n_seq;
  sci=min_en/(sumMFE/n_seq);

  decision_model=get_decision_model(modelDir);

  if (decision_model==NULL){
	nrerror("RNAz 0.1 Error:\nCould not find decision-model. You have to set the RNAZDIR enviroment variable pointing to the model files!\n");
  }

  decValue=999;
  classify(&prob,&decValue,decision_model,id,n_seq,z,sci);
  printf("\n############################  RNAz 0.1  ##############################\n\n");
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
  
  printf("%s",output);

  free(output);
  regression_svm_free();
  svm_destroy_model(decision_model);

  return 0;
}

PRIVATE int read_clustal(FILE *clust, char *AlignedSeqs[], char *names[]) {
   char *line, name[100]={'\0'}, *seq;
   int  n, nn=0, num_seq = 0;
   
   if ((line=get_line(clust)) == NULL) {
     fprintf(stderr, "Empty CLUSTAL file\n"); return 0;
   }

   if (strncmp(line,"CLUSTAL", 7) !=0) {
     fprintf(stderr, "This doesn't look like a CLUSTAL file, sorry\n");
     free(line); return 0;
   }
   free(line);
   line = get_line(clust);

   while (line!=NULL) {
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
       AlignedSeqs[nn] = strdup(seq);
     }
     else {
       if (strcmp(name, names[nn])!=0) {
	 /* name doesn't match */
	 fprintf(stderr,
		 "Sorry, your file is fucked up (inconsitent seq-names)\n");
	 free(line); free(seq);
	 return 0;
       }
       AlignedSeqs[nn] = (char *)
	 xrealloc(AlignedSeqs[nn], strlen(seq)+strlen(AlignedSeqs[nn])+1);
       strcat(AlignedSeqs[nn], seq);
     }
     nn++;
     if (nn>num_seq) num_seq = nn;
     free(seq);
     free(line);
     if (num_seq>=MAX_NUM_NAMES) {
       fprintf(stderr, "Too many sequences in CLUSTAL file");
       return 0;
     }

     line = get_line(clust);
   }
   
   AlignedSeqs[num_seq] = NULL;
   if (num_seq == 0) {
     fprintf(stderr, "No sequences found in CLSUATL file\n");
     return 0;
   }
   n = strlen(AlignedSeqs[0]); 
   for (nn=1; nn<num_seq; nn++) {
     if (strlen(AlignedSeqs[nn])!=n) {
       fprintf(stderr, "Sorry, your file is fucked up.\n"
	       "Unequal lengths!\n\n");
       return 0;
     }
   }
   
   return num_seq;
}
 
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

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

PRIVATE double meanPairID(const char *AS[],int n_seq) {

  int i,j,k,matches,pairs,length;

  matches=0;
  pairs=0;

  length=strlen(AS[0]);
  
  for (i=0;i<n_seq;i++){
	for (j=i+1;j<n_seq;j++){
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

void revAln(char *AS[],int n_seq) {

  int i,j,length;
  char *tmp;
  char letter;
  length=strlen(AS[0]);
  
  for (i=0;i<n_seq;i++){
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

void sliceAln(const char *sourceAln[], char *destAln[], int n_seq, int from, int to){

  int i;
  char *slice;

  for (i=0;i<n_seq;i++){
	slice=(char *) space((unsigned) (to-from+2));
	strncpy(slice,sourceAln[i]+from-1,(to-from+1));
	destAln[i]=slice;
  }
}



PRIVATE void classify(double* prob, double* decValue,struct svm_model* decision_model,double id,int n_seq, double z,double sci){

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


/*-------------------------------------------------------------------------*/

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
