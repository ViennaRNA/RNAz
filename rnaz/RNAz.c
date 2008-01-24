/*********************************************************************                
 *                                                                   *
 *                              RNAz.c                               *
 *                                                                   *
 *	Assess alignments for exceptionally stable and/or conserved      *
 *	secondary structures using RNAalifold/RNAfold and SVMs.          *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                      *
 *                                                                   *
 *	   $Id: RNAz.c,v 1.13 2008-01-24 10:19:14 wash Exp $             *
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
#include "PS_dot.h"
#include "zscore.h"
#include "rnaz_utils.h"
#include "svm.h"
#include "svm_helper.h"
#include "cmdline.h"
#include "strand.h" 

#define IN_RANGE(LOWER,VALUE,UPPER) ((VALUE <= UPPER) && (VALUE >= LOWER))

PRIVATE void usage(void);
PRIVATE void help(void);
PRIVATE void version(void);
PRIVATE void classify(double* prob, double* decValue,
					  struct svm_model* decision_model,
					  double id,int n_seq, double z,double sci);

PRIVATE void gfxPlot(char* fname, char* structure, char* consensus, struct aln *window[],int n_seq);
PRIVATE char **annote(const char *structure, const char *AS[]);


PRIVATE void warning(char* string, double id,int n_seq, double z,double sci,
					 struct aln *AS[]);


enum {FORWARD=1, REVERSE=2};

/********************************************************************
 *                                                                  *
 * main -- main program                                             *
 *                                                                  *
 ********************************************************************/

int main(int argc, char *argv[])
{

  struct svm_model* decision_model; /* SVM classification model */

  /* Command line options */
  int reverse=0;     /* Scan reverse complement */
  int showVersion=0; /* Shows version and exits */
  int showHelp=0;    /* Show short help and exits */
  int from=-1;       /* Scan slice from-to  */
  int to=-1;
  FILE *clust_file=stdin; /* Input file */
  FILE *out=stdout; /* Output file */

  struct aln *AS[MAX_NUM_NAMES];     
  struct aln *window[MAX_NUM_NAMES]; 
  struct aln *windowFwdMem[MAX_NUM_NAMES];

  char *tmpAln[MAX_NUM_NAMES];

  int n_seq;  /* number of input sequences */
  int length; /* length of alignment/window */

  char *structure=NULL;
  char *singleStruc,*gapStruc, *output,*woGapsSeq;
  char strandString[8];
  char warningString[1000];
  char gfxName[100];
  char *string=NULL;
  double singleMFE,sumMFE,singleZ,sumZ,z,sci,id,decValue,prob,comb;
  double min_en, real_en;
  int i,j,k,l,ll,r,countAln;
  int (*readFunction)(FILE *clust,struct aln *alignedSeqs[]);
  char** lines=NULL;
  int directions[3]={FORWARD,0,0};
  int currDirection;
  struct gengetopt_args_info args;

  double meanMFE_fwd=0;
  double consensusMFE_fwd=0;
  double sci_fwd=0;
  double z_fwd=0;
  double gu=0;
  double gu_fwd=0;

  int strandGuess;
  double strandScore,strandProb,strandDec;
  
  if (cmdline_parser (argc, argv, &args) != 0){
	usage();
	exit(EXIT_FAILURE);
  }

  if (args.help_given){
	help();
	exit(EXIT_SUCCESS);
  }

  if (args.version_given){
	version();
	exit(EXIT_SUCCESS);
  }

  if (args.outfile_given){
	out = fopen(args.outfile_arg, "w");
	if (out == NULL){
	  fprintf(stderr, "ERROR: Can't open output file %s\n", args.outfile_arg);
	  exit(1);
	}
  }

  /* Strand prediction implies both strands scored */
  if (args.predict_strand_flag){
	args.both_strands_flag=1;
  }
  
  
  if (args.forward_flag && !args.reverse_flag){
	directions[0]=FORWARD;
	directions[1]=directions[2]=0;
  }
  if (!args.forward_flag && args.reverse_flag){
	directions[0]=REVERSE;
	directions[1]=directions[2]=0;
  }
  if ((args.forward_flag && args.reverse_flag) || args.both_strands_flag){
	directions[0]=FORWARD;
	directions[1]=REVERSE;
  }

  if (args.window_given){
	if (sscanf(args.window_arg,"%d-%d",&from,&to)!=2){
	  nrerror("ERROR: Invalid --window/-w command. "
			  "Use it like '--window 100-200'\n");
	}
	
	printf("from:%d,to:%d\n",from,to);
  }

  
  if (args.inputs_num>=1){
	clust_file = fopen(args.inputs[0], "r"); 
	if (clust_file == NULL){
	  fprintf(stderr, "ERROR: Can't open input file %s\n", args.inputs[0]);
	  exit(1);
	}
  }


 
  /* Global RNA package variables */
  do_backtrack = 1; 
  dangles=2;

  switch(checkFormat(clust_file)){

  case CLUSTAL:
	readFunction=&readClustal;
	break;
  case MAF:
	readFunction=&readMaf;
	break;
  case 0:
	nrerror("ERROR: Unknown alignment file format. Use Clustal W or MAF format.\n");
  }

  decision_model=get_decision_model(NULL);
  regression_svm_init(NULL);
  
  countAln=0;

  while ((n_seq=readFunction(clust_file, AS))!=0){

	if (n_seq ==1){
	  nrerror("ERROR: You need at least two sequences in the alignment.\n");
	}
	
	countAln++;
	
	length = (int) strlen(AS[0]->seq);
	
	/* if a slice is specified by the user */
  
	if ((from!=-1 || to!=-1) && (countAln==1)){
	  
	  if ((from>=to)||(from<=0)||(to>length)){
		nrerror("ERROR: Invalid window range given.\n");
	  }
	  
	  sliceAln((const struct aln**)AS, (struct aln **)window, from, to);
	  length=to-from+1;
	} else { /* take complete alignment */
	  /* window=AS does not work..., deep copy seems not necessary here*/
	  from=1;
	  to=length;
	  sliceAln((const struct aln **)AS, (struct aln **)window, 1, length);
	}

	 /* Convert all Us to Ts for RNAalifold. There is a slight
		 difference in the results. During training we used alignments
		 with Ts, so we use Ts here as well. */

	for (i=0;i<n_seq;i++){
	  j=0;
	  while (window[i]->seq[j]){
		window[i]->seq[j]=toupper(window[i]->seq[j]);
		if (window[i]->seq[j]=='U') window[i]->seq[j]='T';
		++j;
	  }
	}
	
	k=0;
	while ((currDirection=directions[k++])!=0){
	  
	  if (currDirection==REVERSE){
		revAln((struct aln **)window);
		strcpy(strandString,"reverse");
	  } else {
		strcpy(strandString,"forward");
	  }

	  structure = (char *) space((unsigned) length+1);

	  for (i=0;window[i]!=NULL;i++){
		tmpAln[i]=window[i]->seq;
	  }
	  tmpAln[i]=NULL;

	  min_en = alifold(tmpAln, structure);

	  comb=combPerPair(window,structure);
	  
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

		if (window[1]->strand!='?' && !args.window_given){
		  sprintf(output+strlen(output),
				  ">%s %d %d %c %d\n",
				  window[i]->name,window[i]->start,
				  window[i]->length,window[i]->strand,
				  window[i]->fullLength);
		} else {
		  sprintf(output+strlen(output),">%s\n",window[i]->name);
		}

		if (args.show_gaps_flag){

		  gapStruc= (char *) space(sizeof(char)*(strlen(window[i]->seq)+1));

		  l=ll=0;

		  while (window[i]->seq[l]!='\0'){
			if (window[i]->seq[l]!='-'){
			  gapStruc[l]=singleStruc[ll];
			  l++;
			  ll++;
			} else {
			  gapStruc[l]='-';
			  l++;
			}
		  }
		  		  
		  sprintf(output+strlen(output),"%s\n%s ( %6.2f)\n",
				  window[i]->seq,gapStruc,singleMFE);
		  
		  
		} else {
		  sprintf(output+strlen(output),"%s\n%s ( %6.2f)\n",
				  woGapsSeq,singleStruc,singleMFE);
		}
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

	  string = consensusSeq((const struct aln**) window);
	  sprintf(output+strlen(output),
			  ">consensus\n%s\n%s (%6.2f = %6.2f + %6.2f) \n",
			  string, structure, min_en, real_en, min_en-real_en );

	  id=meanPairID((const struct aln**)window);
	  z=sumZ/n_seq;

	  if (sumMFE==0){ 
		/*Set SCI to 0 in the weird case of no structure in single
		  sequences*/
		sci=0;
	  } else {
		sci=min_en/(sumMFE/n_seq);
	  }

	  
	  	  
	  decValue=999;
	  prob=0;
	
	  classify(&prob,&decValue,decision_model,id,n_seq,z,sci);

	  if (args.cutoff_given){
		if (prob<args.cutoff_arg){
		  continue;
		}
	  }

	  
	  strcpy(warningString,"");

	  warning(warningString,id,n_seq,z,sci,(struct aln **)window);
	  

	  fprintf(out,"\n############################  RNAz "PACKAGE_VERSION"  ##############################\n\n");
	  fprintf(out," Sequences: %u\n", n_seq);

	  if (args.window_given){
		fprintf(out," Slice: %u to %u\n",from,to);
	  }
	  fprintf(out," Columns: %u\n",length);
	  fprintf(out," Reading direction: %s\n",strandString);
	  fprintf(out," Mean pairwise identity: %6.2f\n", id);
	  fprintf(out," Mean single sequence MFE: %6.2f\n", sumMFE/n_seq);
	  fprintf(out," Consensus MFE: %6.2f\n",min_en);
	  fprintf(out," Energy contribution: %6.2f\n",real_en);
	  fprintf(out," Covariance contribution: %6.2f\n",min_en-real_en);
	  fprintf(out," Combinations/Pair: %6.2f\n",comb);
	  fprintf(out," Mean z-score: %6.2f\n",z);
	  fprintf(out," Structure conservation index: %6.2f\n",sci);
	  fprintf(out," SVM decision value: %6.2f\n",decValue);
	  fprintf(out," SVM RNA-class probability: %6f\n",prob);
	  if (prob>0.5){
		fprintf(out," Prediction: RNA\n");
	  }
	  else {
		fprintf(out," Prediction: OTHER\n");
	  }
	  fprintf(out,"%s",warningString);
	  
	  fprintf(out,"\n######################################################################\n\n");
	
	  fprintf(out,"%s",output);
	
	  fflush(out);

	  if (args.plot_given){
		if (countAln==1){
		  strcpy(gfxName,".ps");
		} else {
		  sprintf(gfxName,"%i.ps",countAln-1);
		}
		gfxPlot(gfxName, structure, string, window, n_seq);
	  }
 

	  
	 

	  

	  if (currDirection==FORWARD && args.predict_strand_flag){

      sliceAln((const struct aln**)window, (struct aln **)windowFwdMem, 1, length);

      meanMFE_fwd=sumMFE/n_seq;
      consensusMFE_fwd=min_en;
      sci_fwd=sci;
      z_fwd=z;
      
      string = (char*)consensusSeq((const struct aln**) window);
      gu_fwd=(double)frequency_gu_pairs(string, structure);
      
	  }

	  if (currDirection==REVERSE && args.predict_strand_flag){

      string = (char*)consensusSeq((const struct aln**) window);
      gu=(double)frequency_gu_pairs(string, structure);


      // strand(double deltaSCI, double deltaMeanMFE, double deltaConsMFE, double deltaZ, int n_seq, double id, double consFreqGU, 
      //            double *score, double *prob, double *decValue, char *modelDir, int code)

      printf("inhere\n");

      if (strand(sci_fwd-sci, meanMFE_fwd-(sumMFE/n_seq), 
                 consensusMFE_fwd-min_en, z_fwd-z, n_seq, id, gu_fwd+gu,
                 &strandScore, &strandProb, &strandDec, NULL,CODE_SCI_Z_MEANMFE_CONSMFE)){


      /* if (predict_strand(sci_fwd-sci, meanMFE_fwd-(sumMFE/n_seq), */
/*                          consensusMFE_fwd-min_en, z_fwd-z, n_seq, id,  */
/*                          &strandGuess, &strandProb, &strandDec, NULL)){ */

       
        if (strandScore>1){
          fprintf(out, "\n# Strand winner: forward (%.6f)\n",strandDec);
        } else {
          fprintf(out, "\n# Strand winner: reverse (%.6f)\n",strandDec);
        }
      } else {
        fprintf(out, "\n# WARNING: No strand prediction (values out of range)\n");
      }
	  }

    free(structure);
    free(output);

	}
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
/********************************************************************
 *                                                                  *
 * gfxPlot  -- Generates color annotated alignment and RNA          *
 * strucutre as PostScript                                          *  
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * name ... end of filename                                         *
 * structure ... Consensus structure in dot/bracket notation        *
 * consensus ... Consensus sequence                                 *
 * window ... alignment structure                                   *
 * n_seq  ... number of sequences in the alignment                  *
 *                                                                  *
 ********************************************************************/

void gfxPlot(char* name, char* structure, char* consensus, struct aln *window[],int n_seq){ 

  char **A;
  char fname[128];
  int k;
  char **seqs;
  char **names;

  seqs=(char**)space(sizeof(char*)*(n_seq+1));
  names=(char **)space(sizeof(char*)*(n_seq+1));

  for (k=0;k<n_seq;k++){
    seqs[k]=window[k]->seq;
    names[k]=window[k]->name;
  }

  sprintf(fname,"rna%s",name); 

  A = annote(structure, (const char**) seqs);

  (void) PS_rna_plot_a(consensus, structure, fname, A[0], A[1]);

  sprintf(fname,"aln%s",name);

  PS_color_aln(structure, fname,(const char**) seqs,(const char**) names);

  free(A[0]);
  free(A[1]);
  free(seqs);
  free(names);
  
}


PRIVATE char **annote(const char *structure, const char *AS[]) {
  char *ps, *colorps, **A;
  int i, n, s, pairings, maxl;
  short *ptable;
  char * colorMatrix[6][3] = {
    {"0.0 1", "0.0 0.6",  "0.0 0.2"},  /* red    */
    {"0.16 1","0.16 0.6", "0.16 0.2"}, /* ochre  */
    {"0.32 1","0.32 0.6", "0.32 0.2"}, /* turquoise */
    {"0.48 1","0.48 0.6", "0.48 0.2"}, /* green  */
    {"0.65 1","0.65 0.6", "0.65 0.2"}, /* blue   */
    {"0.81 1","0.81 0.6", "0.81 0.2"} /* violet */
  };

  make_pair_matrix();
  n = strlen(AS[0]);
  maxl = 1024;

  A = (char **) space(sizeof(char *)*2);
  ps = (char *) space(maxl);
  colorps = (char *) space(maxl);
  ptable = make_pair_table(structure);
  for (i=1; i<=n; i++) {
    char pps[64], ci='\0', cj='\0';
    int j, type, pfreq[8] = {0,0,0,0,0,0,0,0}, vi=0, vj=0;
    if ((j=ptable[i])<i) continue;
    for (s=0; AS[s]!=NULL; s++) {
      type = pair[encode_char(AS[s][i-1])][encode_char(AS[s][j-1])];
      pfreq[type]++;
      if (type) {
	if (AS[s][i-1] != ci) { ci = AS[s][i-1]; vi++;}
	if (AS[s][j-1] != cj) { cj = AS[s][j-1]; vj++;}
      }
    }
    for (pairings=0,s=1; s<=7; s++) {
      if (pfreq[s]) pairings++;
    }

    if ((maxl - strlen(ps) < 192) || ((maxl - strlen(colorps)) < 64)) {
      maxl *= 2;
      ps = realloc(ps, maxl);
      colorps = realloc(colorps, maxl);
      if ((ps==NULL) || (colorps == NULL))
	  nrerror("out of memory in realloc");
    }

    if (pfreq[0]<=2) {
      snprintf(pps, 64, "%d %d %s colorpair\n",
	       i,j, colorMatrix[pairings-1][pfreq[0]]);
      strcat(colorps, pps);
    }

    if (pfreq[0]>0) {
      snprintf(pps, 64, "%d %d %d gmark\n", i, j, pfreq[0]);
      strcat(ps, pps);
    }
    if (vi>1) {
      snprintf(pps, 64, "%d cmark\n", i);
      strcat(ps, pps);
    }
    if (vj>1) {
      snprintf(pps, 64, "%d cmark\n", j);
      strcat(ps, pps);
    }
  }
  free(ptable);
  A[0]=colorps;
  A[1]=ps;
  return A;
}




/*Hardcoded limits a more sophisticated data-model for meta-model information should be
  considered*/

PRIVATE void warning(char* string, double id,int n_seq, double z,double sci,
					 struct aln *AS[]){

  double GC,A,C;
  int i,j,length,n_A,n_C,n_T,n_G;
  char *seq;

  if (id>100.0) {
	strcpy(string," WARNING: Mean pairwise identity too large.\n");
	string+=strlen(string);
  }

  if (id<52.35) {
	strcpy(string," WARNING: Mean pairwise identity too low.\n");
	string+=strlen(string);
  }

  if (n_seq<2) {
	strcpy(string," WARNING: Too few sequences in alignment.\n");
	string+=strlen(string);
  }
  
  if (n_seq>6) {
	strcpy(string," WARNING: Too many sequences in alignment.\n");
	string+=strlen(string);
  }

  if (!IN_RANGE(-7.87,z,2.76)){
	strcpy(string," WARNING: Mean z-score out of range.\n");
	string+=strlen(string);
  }
  
  if (!IN_RANGE(0,sci,1.23)){
	strcpy(string," WARNING: Structure conservation index out of range.\n");
	string+=strlen(string);
  }

  for (i=0;AS[i]!=NULL;i++){
	seq=AS[i]->seq;

	//printf("SEQ: %s\n",seq);
	
	n_A=n_T=n_G=n_C=0;
  
	for (j=0; j<strlen(seq); j++) {
	  switch (seq[j]) {
	  case 'A': n_A++; break;
	  case 'C': n_C++; break;
	  case 'G': n_G++; break;
	  case 'T': n_T++; break;
	  case 'U': n_T++; break;  
	  }
	}
	
	length=strlen(seq);


	GC=((double)(n_G+n_C)/(double)(n_G+n_C+n_A+n_T));
	A=((double)n_A/(n_A+n_T));
	C=((double)n_C/(n_G+n_C));

	if (length<50){
	  sprintf(string," WARNING: Sequence %d too short.\n",i+1);
	  string+=strlen(string);
	}

	if (length>400){
	  sprintf(string," WARNING: Sequence %d too long.\n",i+1);
	  string+=strlen(string);
	}
	
	if ((!IN_RANGE(0.25,GC,0.75)) ||
		(!IN_RANGE(0.25,GC,0.75)) ||
		(!IN_RANGE(0.25,GC,0.75))){
	  sprintf(string," WARNING: Sequence %d: Base composition out of range.\n",i+1);
	  string+=strlen(string);
	}
  }
}



/********************************************************************
 *                                                                  *
 * usage, help, version - shows information and exits               *
 *                                                                  *
 ********************************************************************/

PRIVATE void usage(void){

  help();

}

PRIVATE void help(void){

  cmdline_parser_print_version ();

  printf("\nUsage: %s [OPTIONS]... [FILES]\n\n", CMDLINE_PARSER_PACKAGE);
  printf("%s\n","  -h, --help              Print help and exit");
  printf("%s\n","  -V, --version           Print version and exit");
  printf("%s\n","  -f, --forward           Score forward strand");
  printf("%s\n","  -r, --reverse           Score reverse strand");
  printf("%s\n","  -b, --both-strands      Score both strands");
  printf("%s\n","  -o, --outfile=FILENAME  Output filename");
  printf("%s\n","  -w, --window=START-STOP Score window from START to STOP");
  printf("%s\n","  -p, --cutoff=FLOAT      Probability cutoff");
  printf("%s\n","  -g, --show-gaps         Display alignment with gap  (default=off)");
  printf("%s\n","  -s, --predict-strand    Use strand predictor  (default=off)");
  printf("%s\n\n","  -x, --plot              Generate graphical output  (default=off)");

}

PRIVATE void version(void){
  printf("RNAz version " PACKAGE_VERSION "\n");
  exit(EXIT_SUCCESS);
}







