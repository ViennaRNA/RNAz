/*********************************************************************                
 *                                                                   *
 *                              RNAz.c                               *
 *                                                                   *
 *	Assess alignments for exceptionally stable and/or conserved      *
 *	secondary structures using RNAalifold/RNAfold and SVMs.          *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                      *
 *                                                                   *
 *	   $Id: RNAz.c,v 1.7 2006-01-09 18:32:52 wash Exp $              *
 *                                                                   *
 *********************************************************************/
#include "config.h"
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

enum alnFormat {UNKNOWN=0, CLUSTAL=1, MAF=2};

struct aln {
  char *name;
  char *seq;
};


PRIVATE void usage(void);
PRIVATE void help(void);
PRIVATE void version(void);

PRIVATE int read_clustal(FILE *clust,
						 struct aln *alignedSeqs[]);

PRIVATE int read_maf(FILE *clust,
						 struct aln *alignedSeqs[]);


PRIVATE char *consensus(const struct aln *AS[]);

PRIVATE double meanPairID(const struct aln *AS[]);

PRIVATE void revAln(struct aln *AS[]);

PRIVATE void sliceAln(const struct aln *sourceAln[], struct aln *destAln[],
					  int from, int to);

PRIVATE void freeAln(struct aln *AS[]);

PRIVATE void classify(double* prob, double* decValue,
					  struct svm_model* decision_model,
					  double id,int n_seq, double z,double sci);

PRIVATE struct aln* createAlnEntry(char* name, char* seq); 
PRIVATE void freeAlnEntry(struct aln* entry);
PRIVATE void printAln(const struct aln* AS[]);

PRIVATE int checkFormat(FILE *file);

PRIVATE char** splitFields(char* string);
PRIVATE void freeFields(char** fields);


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
  char  *outFile=NULL;
  FILE *out=stdout; /* Output file */

  struct aln *AS[MAX_NUM_NAMES];     
  struct aln *window[MAX_NUM_NAMES]; 
  char *tmpAln[MAX_NUM_NAMES];

  int n_seq;  /* number of input sequences */
  int length; /* length of alignment/window */

  char *structure=NULL;
  char *singleStruc, *output,*woGapsSeq;
  char strand[8];
  char *string=NULL;
  double singleMFE,sumMFE,singleZ,sumZ,z,sci,id,decValue,prob;
  double min_en, real_en;
  int i,j,r,countAln;
  int (*readFunction)(FILE *clust,struct aln *alignedSeqs[]);

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
	   case 'o':
		 outFile = argv[++i];
		 if (outFile[0]=='-'){
		   usage();
		 }
		 out = fopen(outFile, "w");
		 if (out == NULL){
		   nrerror("Could not open output file");
		   exit(1);
		 }
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
  
  switch(checkFormat(clust_file)){

  case CLUSTAL:
	readFunction=&read_clustal;
	break;
  case MAF:
	readFunction=&read_maf;
	break;
  }

  if (showVersion==1) version();
  if (showHelp==1) help();
  
  modelDir=getenv("RNAZDIR");


  if (modelDir==NULL){
	nrerror("ERROR: Could not find models. You have to set the RNAZDIR"
			" enviroment variable pointing to the model files!\n");
  }
  

  regression_svm_init(modelDir);
  
  decision_model=get_decision_model(modelDir);
  

  if (decision_model==NULL){	
	nrerror("ERROR: Could not find decision-model. You have to set "			
			"the RNAZDIR enviroment variable pointing to the model files!\n");	
  }

  
  countAln=0;
  while ((n_seq=readFunction(clust_file, AS))!=0){
	countAln++;
	length = (int) strlen(AS[0]->seq);
	
	/* if a slice is specified by the user */
  
	if ((from!=-1 || to!=-1) && (countAln==1)){
	  if (((from!=-1) && (to==-1)) || ((from==-1) && (to!=-1))){
		nrerror("ERROR: You have to set both -f and -t parameters"
				" to score a slice of the alignment.\n");
	  }
	  if ((from>=to)||(from<=0)||(to>length)){
		nrerror("ERROR: -f and/or -t parameters out of range"
				" (no reasonable slice specified)\n");
	  }
	  sliceAln((const struct aln**)AS, (struct aln **)window, from, to);
	  length=to-from+1;
	} else { /* take complete alignment */
	  /* window=AS does not work..., deep copy seems not necessary here*/
	  from=1;
	  to=length;
	  sliceAln((const struct aln **)AS, (struct aln **)window, 1, length);
	}
  
	if (reverse==1){
	  revAln((struct aln **)window);
	  strcpy(strand,"reverse");
	} else {
	  strcpy(strand,"forward");
	}

	//printAln((const struct aln**)window);
	
	structure = (char *) space((unsigned) length+1);

	for (i=0;window[i]!=NULL;i++){
	  tmpAln[i]=window[i]->seq;
	}
	tmpAln[i]=NULL;



	/* Convert all Us to Ts for RNAalifold. There is a slight
	   difference in the results. During training we used alignments
	   with Ts, so we use Ts here as well. */

	for (i=0;i<n_seq;i++){
	  j=0;
	  while (window[i]->seq[j]){
		window[i]->seq[j]=toupper(window[i]->seq[j]);
		if (window[i]->seq[j]=='U') window[i]->seq[j]='T'; ++j;
	  }
	}
	
	min_en = alifold(tmpAln, structure);

	sumZ=0;
	sumMFE=0;

	output=(char *)space(sizeof(char)*(length+16)*(n_seq+1)*3);
  
	for (i=0;i<n_seq;i++){
	  singleStruc = space(strlen(window[i]->seq)+1);
	  woGapsSeq = space(strlen(window[i]->seq)+1);
	  j=0;
	  while (window[i]->seq[j]){
		/* Convert all Ts to Us for RNAfold. There is a difference
		   between the results. With U in the function call, we get
		   the results as RNAfold gives on the command line. Since
		   this variant was also used during training, we use it here
		   as well. */
		if (window[i]->seq[j]=='T') window[i]->seq[j]='U';
		if (window[i]->seq[j]!='-'){
		  woGapsSeq[strlen(woGapsSeq)]=window[i]->seq[j];
		  woGapsSeq[strlen(woGapsSeq)]='\0';
		}
		++j;
	  }
	  
	  singleMFE = fold(woGapsSeq, singleStruc);
	  singleZ=mfe_zscore(woGapsSeq,singleMFE);
	  
	  sumZ+=singleZ;
	  sumMFE+=singleMFE;
	  sprintf(output+strlen(output),">%s\n%s\n%s ( %6.2f)\n",
			  window[i]->name,woGapsSeq,singleStruc,singleMFE);
	  free(woGapsSeq);
	  free(singleStruc);
	}

	{
	  int i; double s=0;
	  extern int eos_debug;
 	  eos_debug=-1; /* shut off warnings about nonstandard pairs */
	  for (i=0; window[i]!=NULL; i++) 
		s += energy_of_struct(window[i]->seq, structure);
	  real_en = s/i;
	}

	string = consensus((const struct aln**) window);
	sprintf(output+strlen(output),
			">consensus\n%s\n%s (%6.2f = %6.2f + %6.2f) \n",
			string, structure, min_en, real_en, min_en-real_en );

	id=meanPairID((const struct aln**)window);
	z=sumZ/n_seq;
	sci=min_en/(sumMFE/n_seq);

	decValue=999;
	prob=0;
	
	classify(&prob,&decValue,decision_model,id,n_seq,z,sci);

	
	fprintf(out,"\n###########################  RNAz "PACKAGE_VERSION"  #############################\n\n");
	fprintf(out," Sequences: %u\n", n_seq);
	fprintf(out," Slice: %u to %u\n",from,to);
	fprintf(out," Columns: %u\n",length);
	fprintf(out," Strand: %s\n",strand);
	fprintf(out," Mean pairwise identity: %6.2f\n", id);
	fprintf(out," Mean single sequence MFE: %6.2f\n", sumMFE/n_seq);
	fprintf(out," Consensus MFE: %6.2f\n",min_en);
	fprintf(out," Energy contribution: %6.2f\n",real_en);
	fprintf(out," Covariance contribution: %6.2f\n",min_en-real_en);
	fprintf(out," Mean z-score: %6.2f\n",z);
	fprintf(out," Structure conservation index: %6.2f\n",sci);
	fprintf(out," SVM decision value: %6.2f\n",decValue);
	fprintf(out," SVM RNA-class probability: %6f\n",prob);
	if (prob>0.5){
	  fprintf(out," Prediction: RNA\n");
	}
	else {
	  fprintf(out," Prediction: no RNA\n");
	}
	
	fprintf(out,"\n######################################################################\n\n");
	
	fprintf(out,"%s//\n\n",output);
	
	fflush(out);

	free(structure);
	free(output);
	freeAln((struct aln **)AS);
	freeAln((struct aln **)window);
	

  }

  if (countAln==0){
	nrerror("ERROR: Empty alignment file\n");
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


PRIVATE int read_clustal(FILE *clust, struct aln *alignedSeqs[]) {

  char *line, name[100]={'\0'}, *seq;
  int  n, nn=0, num_seq = 0;

  if (feof(clust)){
	return 0;
  }
  
  line = get_line(clust);

  while (line!=NULL) {

	if (strncmp(line,"CLUSTAL", 7)==0) {
	  break;
	}
	
	if (((n=strlen(line))<4) || isspace((int)line[0])) {
	  /* skip non-sequence line */
	  free(line); line = get_line(clust);
	  nn=0; /* reset seqence number */
	  continue;
	} 
     
	seq = (char *) space( (n+1)*sizeof(char) );
	sscanf(line,"%99s %s", name, seq);
	if (nn == num_seq) { /* first time */
	  //names[nn] = strdup(name);
	  //alignedSeqs[nn] = strdup(seq);

	  alignedSeqs[nn]=createAlnEntry(strdup(name),strdup(seq));
	  
	  
	}
	else {
	  if (strcmp(name, alignedSeqs[nn]->name)!=0) {
		/* name doesn't match */
		free(line); free(seq);
		nrerror("ERROR: Inconsistent sequence names in CLUSTAL file");
		return 0;
	  }
	  alignedSeqs[nn]->seq = (char *)
		xrealloc(alignedSeqs[nn]->seq, strlen(seq)+strlen(alignedSeqs[nn]->seq)+1);
	  strcat(alignedSeqs[nn]->seq, seq);
	}
	nn++;
	if (nn>num_seq) num_seq = nn;
	free(seq);
	free(line);
	if (num_seq>=MAX_NUM_NAMES) {
	  nrerror("ERROR: Too many sequences in CLUSTAL file");
	  return 0;
	}
	line = get_line(clust);
  }

  alignedSeqs[num_seq] = NULL;

  if (num_seq == 0) {
	return 0;
  }
  n = strlen(alignedSeqs[0]->seq); 
  for (nn=1; nn<num_seq; nn++) {
	if (strlen(alignedSeqs[nn]->seq)!=n) {
	  fprintf(stderr, "ERROR: Sequences are of unequal length.\n");
	  return 0;
	}
  }
  return num_seq;
}

/********************************************************************
 *                                                                  *
 * read_maf -- read MAF formatted file                              *
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


PRIVATE int read_maf(FILE *clust, struct aln *alignedSeqs[]) {

  char *line, name[100]={'\0'}, *seq;
  int num_seq = 0;
  char** fields;
  int n,nn;

  if (feof(clust)){
	return 0;
  }
  
  while ((line=get_line(clust))!=NULL) {

	fields=splitFields(line);

	/* Skip empty (=only whitespace) lines */
	if (fields==NULL){
	  free(line);
	  continue;
	}

	/* Skip comment (#) lines */
	if (fields[0][0]=='#'){
	  free(line);
	  freeFields(fields);
	  continue;
	}

	if (fields[0][0]=='s' && fields[0][1]=='\0'){

	  n=0;
	  while (fields[n]!=NULL){
		n++;
	  }
	  if (n!=7){
		nrerror("ERROR: Invalid MAF format (number of fields in 's' line not correct)");
	  }
	  
	  alignedSeqs[num_seq++]=createAlnEntry(strdup(fields[1]),strdup(fields[6]));
	  free(line);
	  freeFields(fields);
	  continue;
	}

	if (fields[0][0]=='a' && fields[0][1]=='\0'){
	  free(line);
	  freeFields(fields);
	  break;
	}
  }

  alignedSeqs[num_seq] = NULL;

  n = strlen(alignedSeqs[0]->seq); 
  for (nn=1; nn<num_seq; nn++) {
	if (strlen(alignedSeqs[nn]->seq)!=n) {
	  nrerror("ERROR: Sequences are of unequal length.");
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

PRIVATE char *consensus(const struct aln *AS[]) {
  char *string;
  int i,n;
  n = strlen(AS[0]->seq);
  string = (char *) space((n+1)*sizeof(char));
  for (i=0; i<n; i++) {
    int s,c,fm, freq[8] = {0,0,0,0,0,0,0,0};
    for (s=0; AS[s]!=NULL; s++) 
      freq[encode_char(AS[s]->seq[i])]++;
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


PRIVATE double meanPairID(const struct aln *AS[]) {

  int i,j,k,matches,pairs,length;

  matches=0;
  pairs=0;

  length=strlen(AS[0]->seq);
  
  for (i=0;AS[i]!=NULL;i++){
	for (j=i+1;AS[j]!=NULL;j++){
	  for (k=0;k<length;k++){
		if ((AS[i]->seq[k]!='-') || (AS[j]->seq[k]!='-')){
		  if (AS[i]->seq[k]==AS[j]->seq[k]){
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

void revAln(struct aln *AS[]) {

  int i,j,length;
  char *tmp;
  char letter;
  length=strlen(AS[0]->seq);
  
  for (i=0;AS[i]!=NULL;i++){
	tmp = (char *) space((unsigned) length+1);
	for (j=length-1;j>=0;j--){
	  letter=AS[i]->seq[j];
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
	strcpy(AS[i]->seq,tmp);
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

void sliceAln(const struct aln *sourceAln[], struct aln *destAln[],
			  int from, int to){

  int i;
  char *slice;
 
  for (i=0;sourceAln[i]!=NULL;i++){
	slice=(char *) space((unsigned) (to-from+2));
	strncpy(slice,(sourceAln[i]->seq)+from-1,(to-from+1));
	destAln[i]=createAlnEntry(strdup(sourceAln[i]->name),slice);
  }
  destAln[i]=NULL;
}

/********************************************************************
 *                                                                  *
 * freeAln -- Frees memory of alignment array                       *
 *                                                                  *
 ********************************************************************/

PRIVATE void freeAln(struct aln *AS[]){

  int i;
  for (i=0;AS[i]!=NULL;i++){
	freeAlnEntry(AS[i]);
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


  double* value;
  value=(double*)space(sizeof(double)*2);

  node[0].index = 1; node[0].value = z;
  node[1].index = 2; node[1].value = sci;
  node[2].index = 3; node[2].value = id;
  node[3].index = 4; node[3].value = n_seq;
  node[4].index =-1;

  scale_decision_node((struct svm_node*)&node);

  svm_predict_values(decision_model,node,value);
  *decValue=value[0];

  svm_predict_probability(decision_model,node,value);
  *prob=value[0];
  
  
}

PRIVATE struct aln* createAlnEntry(char* name, char* seq){

  struct aln* entry;

  entry=(struct aln*)space(sizeof(struct aln));

  entry->name=name;
  entry->seq=seq;

  return entry;
}

PRIVATE void freeAlnEntry(struct aln* entry){

  free(entry->name);
  free(entry->seq);
  free(entry);
 
}


PRIVATE void printAln(const struct aln* AS[]){
  int i;
  for (i=0;AS[i]!=NULL;i++){
	printf("%s %s\n",AS[i]->name,AS[i]->seq);
  }
}

PRIVATE int checkFormat(FILE *file){

  char *line; 
  char** fields;
  
  while ((line=get_line(file)) != NULL){

	fields=splitFields(line);

	/* Skip empty (=only whitespace) and comments */

	if (fields==NULL){
	  free(line);
	  freeFields(fields);
	  continue;
	}
	
	if (fields[0][0]=='#'){
	  free(line);
	  freeFields(fields);
	  continue;
	}

	/* Identify "CLUSTAL" header => CLUSTAL*/
	if (strcmp(fields[0],"CLUSTAL")==0){
	  free(line);
	  freeFields(fields);
	  return(CLUSTAL);
	}

	/* Identitfy "a" header => MAF */
	if (fields[0][0]=='a' && fields[0][1]=='\0'){
	  free(line);
	  freeFields(fields);
	  return(MAF);
	}

	/* Unknown format if the first non-empty, non-command data in the
	   file is non of the above */
	free(line);
	freeFields(fields);
	return(0);
  }

  free(line);
  nrerror("ERROR: Empty alignment file\n");
   
}



PRIVATE char** splitFields(char* string){

  char c;
  char* currField;
  char** output=NULL;
  int* seps;
  int nSep;
  int nField=0;
  int i=0;

  if (strlen(string)==0 || string==NULL){
	return NULL;
  }

  /* First find all characters which are whitespaces and store the
	 positions in the array seps */
    
  seps=(int *)space(sizeof(int));
  seps[0]=-1;
  nSep=1;
  
  while ((c=string[i])!='\0' && (c!='\n')){
	if (isspace(c)){
	  seps=(int*)xrealloc(seps,sizeof(int)*(nSep+1));
	  seps[nSep++]=i;
	}
	i++;
  }

  seps=(int*)xrealloc(seps,sizeof(int)*(nSep+1));
  seps[nSep]=strlen(string);


  /* Then go through all intervals in between of two whitespaces (or
	 end or start of string) and store the fields in the array
	 "output"; if there are two adjacent whitespaces this is ignored
	 resulting in a behaviour like "split /\s+/" in perl */
  
  for (i=0;i<nSep;i++){

	int start=seps[i];
	int stop=seps[i+1];
	int length=(stop-start);
	int notSpace,j;

	
	currField=(char *)space(sizeof(char)*(length+1));
	strncpy(currField,string+start+1,length-1);
	currField[length]='\0';

	/* check if field is not only whitespace */
	notSpace=0;
	j=0;
	while (c=currField[j]!='\0'){
	  if (!isspace(c)){
		notSpace=1;
		break;
	  }
	}

	if (notSpace){
	  output=(char**)xrealloc(output,sizeof(char**)*(nField+1));
	  output[nField++]=currField;
	  currField=NULL;
	} else {
	  free(currField);
	  currField=NULL;
	}

	//printf("%s|\n",output[nField-1]);
  }

  if (nField==0){
	return NULL;
  }

  
  output=(char**)xrealloc(output,sizeof(char**)*(nField+1));
  output[nField]=NULL;
  
  free(seps);
  return output;
  
}

PRIVATE void freeFields(char** fields){

  int i=0;
  while (fields[i]!=NULL){
	free(fields[i++]);
  }
  free(fields);
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
		  "-o      Write output to file (default STDOUT)\n"
		  "-v      Show version information\n"
		  "-h      Show this help message\n"
		  );

}

PRIVATE void version(void){
  printf("RNAz version " PACKAGE_VERSION ", September 2004\n");
  exit(EXIT_SUCCESS);
}

