/*********************************************************************                
 *                                                                   *
 *                              RNAz.c                               *
 *                                                                   *
 *	Assess alignments for exceptionally stable and/or conserved      *
 *	secondary structures using RNAalifold/RNAfold and SVMs.          *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                      *
 *                                                                   *
 *	   $Id: RNAz.c,v 1.12 2006/10/12 16:58:18 wash Exp $             *
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
		      double id,int n_seq, double z,double sci,
		      double entropy, int decision_model_type);
PRIVATE void warning(char* string, double id, int n_seq, 
		     double z, double sci, double entropy,
		     struct aln *AS[], int decision_model_type);


enum {FORWARD=1, REVERSE=2};

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
  FILE *out=stdout; /* Output file */

  struct aln *AS[MAX_NUM_NAMES];     
  struct aln *window[MAX_NUM_NAMES]; 
  char *tmpAln[MAX_NUM_NAMES];

  int n_seq;  /* number of input sequences */
  int length; /* length of alignment/window */
  int z_score_type;
  int decision_model_type;

  char *structure=NULL;
  char *singleStruc,*gapStruc, *output,*woGapsSeq;
  char strand[8];
  char warningString[2000];
  char warningString_regression[2000];
  char *string=NULL;
  double singleMFE,sumMFE,singleZ,sumZ,z,sci,id,decValue,prob,comb,entropy,GC;
  double min_en, real_en;
  int i,j,k,l,ll,r,countAln,nonGaps,singleGC;
  int (*readFunction)(FILE *clust,struct aln *alignedSeqs[]);
  char** lines=NULL;
  int directions[3]={FORWARD,0,0};
  int currDirection;
  struct gengetopt_args_info args;

  double meanMFE_fwd=0;
  double consensusMFE_fwd=0;
  double sci_fwd=0;
  double z_fwd=0;
  int strandGuess;
  int avoid_shuffle=0;
  double strandProb,strandDec;
  
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
	readFunction=&read_clustal;
	break;
  case MAF:
	readFunction=&read_maf;
	break;
  case 0:
	nrerror("ERROR: Unknown alignment file format. Use Clustal W or MAF format.\n");
  }

  /* Set z-score type (mono/dinucleotide) here */
  z_score_type = 0;
  if (args.dinucleotide_given) z_score_type = 2;

  /* now let's decide which decision model to take */
  /* decision_model_type = 1 for normal model used in RNAz 1.0 */
  /* decision_model_type = 2 for normal model using dinucelotide background */
  /* decision_model_type = 3 for structural model using dinucelotide background */
  decision_model_type = 1;
  if (args.dinucleotide_given) decision_model_type = 2;
  if (args.dinucleotide_given && args.locarnate_given) decision_model_type = 3;
  if ((!args.dinucleotide_given) && args.locarnate_given){
    z_score_type=2;
    //nrerror("ERROR: Structural decision model only trained with dinucleotide background model.\n");
  }


  if (args.no_shuffle_given) avoid_shuffle = 1;

  decision_model=get_decision_model(NULL, decision_model_type);

  /* Initialize Regression Models for mononucleotide */
  /* Not needed if we score with dinucleotides */
  if (z_score_type == 0) regression_svm_init();
  
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
		strcpy(strand,"reverse");
	  } else {
		strcpy(strand,"forward");
	  }

	  structure = (char *) space((unsigned) length+1);

	  for (i=0;window[i]!=NULL;i++){
		tmpAln[i]=window[i]->seq;
	  }
	  tmpAln[i]=NULL;

	  min_en = alifold(tmpAln, structure);
	  free_alifold_arrays();

	  comb=combPerPair(window,structure);
	  
	  sumZ=0.0;
	  sumMFE=0.0;
	  GC=0.0;

	  output=(char *)space(sizeof(char)*(length+160)*(n_seq+1)*3);

	  strcpy(warningString,"");
	  strcpy(warningString_regression,"");

	  for (i=0;i<n_seq;i++){
		singleStruc = space(strlen(window[i]->seq)+1);
		woGapsSeq = space(strlen(window[i]->seq)+1);
		j=0;
		nonGaps=0;
		singleGC=0;
		while (window[i]->seq[j]){
		  /* Convert all Ts to Us for RNAfold. There is a difference
		     between the results. With U in the function call, we get
		     the results as RNAfold gives on the command line. Since
		     this variant was also used during training, we use it here
		     as well. */
		  if (window[i]->seq[j]=='T') window[i]->seq[j]='U';
		  if (window[i]->seq[j]=='C') singleGC++;
		  if (window[i]->seq[j]=='G') singleGC++;
		  if (window[i]->seq[j]!='-'){
		    nonGaps++;
		    woGapsSeq[strlen(woGapsSeq)]=window[i]->seq[j];
		    woGapsSeq[strlen(woGapsSeq)]='\0';
		  }
		  ++j;
		}
		
		/* z-score is calculated here! */
		singleMFE = fold(woGapsSeq, singleStruc);
		free_arrays();
		/* z-score type may be overwritten. If it is out of training
		   bounds, we switch to shuffling if allowed (avoid_shuffle). */
		int z_score_type_orig = z_score_type;


		singleZ=mfe_zscore(woGapsSeq,singleMFE, &z_score_type, avoid_shuffle, warningString_regression);

		GC+=(double) singleGC/nonGaps;
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
		  char ch;
		  ch = 'R';
		  if (z_score_type == 1 || z_score_type == 3) ch = 'S';
		  		  
		  sprintf(output+strlen(output),"%s\n%s ( %6.2f, z-score = %6.2f, %c)\n",
			  window[i]->seq,gapStruc,singleMFE,singleZ,ch);
		  z_score_type = z_score_type_orig;
		  
		} else {
		  char ch;
		  ch = 'R';
		  if (z_score_type == 1 || z_score_type == 3) ch = 'S';

		  sprintf(output+strlen(output),"%s\n%s ( %6.2f, z-score = %6.2f, %c)\n",
			  woGapsSeq,singleStruc,singleMFE,singleZ,ch);
		  z_score_type = z_score_type_orig;
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

	  string = consensus((const struct aln**) window);
	  sprintf(output+strlen(output),
			  ">consensus\n%s\n%s (%6.2f = %6.2f + %6.2f) \n",
			  string, structure, min_en, real_en, min_en-real_en );
	  free(string);

	  id=meanPairID((const struct aln**)window);
	  entropy=NormShannonEntropy((const struct aln**)window);
	  z=sumZ/n_seq;
	  GC=(double)GC/n_seq;

	  if (sumMFE==0){ 
		/*Set SCI to 0 in the weird case of no structure in single
		  sequences*/
		sci=0;
	  } else {
		sci=min_en/(sumMFE/n_seq);
	  }

	  
	  	  
	  decValue=999;
	  prob=0;
	
	  classify(&prob,&decValue,decision_model,id,n_seq,z,sci,entropy,decision_model_type);

	  if (args.cutoff_given){
		if (prob<args.cutoff_arg){
		  continue;
		}
	  }

	  warning(warningString,id,n_seq,z,sci,entropy,(struct aln **)window,decision_model_type);
	  

 	  fprintf(out,"\n############################  RNAz "PACKAGE_VERSION"  ##############################\n\n"); 
 	  fprintf(out," Sequences: %u\n", n_seq); 

 	  if (args.window_given){ 
 		fprintf(out," Slice: %u to %u\n",from,to); 
 	  } 
 	  fprintf(out," Columns: %u\n",length);
 	  fprintf(out," Reading direction: %s\n",strand); 
	  fprintf(out," Mean pairwise identity: %6.2f\n", id);
	  fprintf(out," Shannon entropy: %2.5f\n", entropy);
	  fprintf(out," G+C content: %2.5f\n", GC);
 	  fprintf(out," Mean single sequence MFE: %6.2f\n", sumMFE/n_seq); 
 	  fprintf(out," Consensus MFE: %6.2f\n",min_en); 
 	  fprintf(out," Energy contribution: %6.2f\n",real_en); 
 	  fprintf(out," Covariance contribution: %6.2f\n",min_en-real_en); 
 	  fprintf(out," Combinations/Pair: %6.2f\n",comb); 
	  fprintf(out," Mean z-score: %6.2f\n",z);
	  fprintf(out," Structure conservation index: %6.2f\n",sci);
	  if (decision_model_type == 1) {
	    fprintf(out," Background model: mononucleotide\n");
	    fprintf(out," Decision model: sequence based alignment quality\n");
	  }
	  if (decision_model_type == 2) {
	    fprintf(out," Background model: dinucleotide\n");
	    fprintf(out," Decision model: sequence based alignment quality\n");
	  }
	  if (decision_model_type == 3) {
	    fprintf(out," Background model: dinucleotide\n");
	    fprintf(out," Decision model: structural RNA alignment quality\n");
	  }
 	  fprintf(out," SVM decision value: %6.2f\n",decValue); 
 	  fprintf(out," SVM RNA-class probability: %6f\n",prob); 
 	  if (prob>0.5){ 
 		fprintf(out," Prediction: RNA\n"); 
 	  } 
 	  else { 
 		fprintf(out," Prediction: OTHER\n"); 
 	  } 

	  fprintf(out,"%s",warningString_regression);

	  fprintf(out,"%s",warningString);
	  
 	  fprintf(out,"\n######################################################################\n\n"); 
	
 	  fprintf(out,"%s",output); 
	
	  fflush(out);

	  free(structure);
	  free(output);

	  if (currDirection==FORWARD && args.predict_strand_flag){
		meanMFE_fwd=sumMFE/n_seq;
		consensusMFE_fwd=min_en;
		sci_fwd=sci;
		z_fwd=z;
	  }

	  if (currDirection==REVERSE && args.predict_strand_flag){

		if (predict_strand(sci_fwd-sci, meanMFE_fwd-(sumMFE/n_seq),
						   consensusMFE_fwd-min_en, z_fwd-z, n_seq, id, 
						   &strandGuess, &strandProb, &strandDec, NULL)){
		  if (strandGuess==1){
			fprintf(out, "\n# Strand winner: forward (%.2f)\n",strandProb);
		  } else {
			fprintf(out, "\n# Strand winner: reverse (%.2f)\n",1-strandProb);
		  }
		} else {
		  fprintf(out, "\n# WARNING: No strand prediction (values out of range)\n");
		}
	  }
	}
	freeAln((struct aln **)AS);
	freeAln((struct aln **)window);
	
  }
  if (args.inputs_num>=1){
    fclose(clust_file);
  }
  cmdline_parser_free (&args);

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
 * entropy ... normalized Shannon entropy                           *
 * decision_model_type ... mono/dinuclotide or structural           *
 *                                                                  *
 ********************************************************************/


PRIVATE void classify(double* prob, double* decValue,
		      struct svm_model* decision_model,
		      double id,int n_seq, double z,double sci,
		      double entropy, int decision_model_type){

  FILE *out=stdout; /* Output file */

  /************************************/
  /* normal model as used in RNAz 1.0 */
  /************************************/
  if (decision_model_type == 1) {
    struct svm_node node[5];
    
    double* value;
    value=(double*)space(sizeof(double)*2);
    
    node[0].index = 1; node[0].value = z;
    node[1].index = 2; node[1].value = sci;
    node[2].index = 3; node[2].value = id;
    node[3].index = 4; node[3].value = n_seq;
    node[4].index =-1;
    
    scale_decision_node((struct svm_node*)&node,decision_model_type);
    
    svm_predict_values(decision_model,node,value);
    *decValue=value[0];
    
    svm_predict_probability(decision_model,node,value);
    *prob=value[0];

    free(value);
  }

  /****************************************************/
  /* dinucleotide model for sequence based alignments */
  /****************************************************/
  if (decision_model_type == 2) {

    /* For training we used the z-score and the SCI rounded to tow 
       decimal places. To make results comparable we do it also here. */

    char tmp[10];
    sprintf(tmp, "%.2f", (double) z);
    double z_tmp = (double) atof(tmp);

    sprintf(tmp, "%.2f", (double) sci);
    double sci_tmp = (double) atof(tmp);

    /* In some rare cases it might happen that z-score and SCI are
       out of the training range. In these cases they are just set to the 
       maximum or minimum. */
 
    if (z_tmp > 2.01) z_tmp = 2.01;
    if (z_tmp < -8.15) z_tmp = -8.15;
    if (sci_tmp > 1.29) sci_tmp = 1.29;

    /* construct SVM node and scale parameters */

    struct svm_node node[4];
    double* value;
    value=(double*)space(sizeof(double)*2);
    
    node[0].index = 1; node[0].value = z_tmp;
    node[1].index = 2; node[1].value = sci_tmp;
    node[2].index = 3; node[2].value = entropy;
    node[3].index =-1; 
    
    scale_decision_node((struct svm_node*)&node,decision_model_type);

    /* For training we used scaled variables rounded to five
       decimal places. To make results comparable we do it also here. */

    sprintf(tmp, "%.5f", (double) node[0].value);
    node[0].value = (double) atof(tmp);
    sprintf(tmp, "%.5f", (double) node[1].value);
    node[1].value = (double) atof(tmp);
    sprintf(tmp, "%.5f", (double) node[2].value);
    node[2].value = (double) atof(tmp);

    /* Now predict decision value and probability */

    svm_predict_values(decision_model,node,value);
    *decValue=value[0];
    
    svm_predict_probability(decision_model,node,value);
    *prob=value[0];

    free(value);
  }


  /***********************************************************************/
  /* dinucleotide model for structural alignments generated by locarnate */
  /***********************************************************************/
  if (decision_model_type == 3) {

    /* For training we used the z-score and the SCI rounded to tow 
       decimal places. To make results comparable we do it also here. */

    char tmp[10];
    sprintf(tmp, "%.2f", (double) z);
    double z_tmp = (double) atof(tmp);

    sprintf(tmp, "%.2f", (double) sci);
    double sci_tmp = (double) atof(tmp);

    /* In some rare cases it might happen that z-score and SCI are
       out of the training range. In these cases they are just set to the 
       maximum or minimum. */
 
    if (z_tmp > 2.01) z_tmp = 2.01;
    if (z_tmp < -8.13) z_tmp = -8.13;
    if (sci_tmp > 1.31) sci_tmp = 1.31;

    /*fprintf(out," z: %f, sci: %f, entropy:%f\n",z_tmp,sci_tmp,entropy);*/

    /* construct SVM node and scale parameters */

    struct svm_node node[4];
    double* value;
    value=(double*)space(sizeof(double)*2);
    
    node[0].index = 1; node[0].value = z_tmp;
    node[1].index = 2; node[1].value = sci_tmp;
    node[2].index = 3; node[2].value = entropy;
    node[3].index =-1; 
    
    scale_decision_node((struct svm_node*)&node,decision_model_type);

    /* For training we used scaled variables rounded to five
       decimal places. To make results comparable we do it also here. */

    sprintf(tmp, "%.5f", (double) node[0].value);
    node[0].value = (double) atof(tmp);
    sprintf(tmp, "%.5f", (double) node[1].value);
    node[1].value = (double) atof(tmp);
    sprintf(tmp, "%.5f", (double) node[2].value);
    node[2].value = (double) atof(tmp);

    /*fprintf(out," z: %f, sci: %f, entropy:%f\n",node[0].value,node[1].value,node[2].value);*/

    /* Now predict decision value and probability */

    svm_predict_values(decision_model,node,value);
    *decValue=value[0];
    
    svm_predict_probability(decision_model,node,value);
    *prob=value[0];

    free(value);
  }
  
  
}

/*Hardcoded limits a more sophisticated data-model for meta-model information should be
  considered*/

PRIVATE void warning(char* string, double id, int n_seq, 
		     double z, double sci, double entropy,
		     struct aln *AS[], int decision_model_type){

  /* Now we throw warnings fors the old RNAz 1.0 */
  
  if (decision_model_type == 1) {
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
	  (!IN_RANGE(0.25,A,0.75)) ||
	  (!IN_RANGE(0.25,C,0.75))){
	sprintf(string," WARNING: Sequence %d: Base composition out of range.\n",i+1);
	string+=strlen(string);
      }
    }
  }

  /* Now we throw warnings for the new RNAz 2.0. There are no limits anymore for 
     z-score regression, so we do not have to care about length or base composition.
     We can also ignore number of sequences and MPI, just check the entropy. */
  
  if (decision_model_type == 2) {
    
    if (!IN_RANGE(-8.15,z,2.01)){
      strcpy(string," WARNING: Mean z-score out of range.\n");
      string+=strlen(string);
    }
    
    if (!IN_RANGE(0,sci,1.29)){
      strcpy(string," WARNING: Structure conservation index out of range.\n");
      string+=strlen(string);
    }

    if (!IN_RANGE(0,entropy,1.28718)){
      strcpy(string," WARNING: Normalized Shannon entropy out of range. Your alignment has too much sequence variation.\n");
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
  printf("%s\n","  -d, --dinucleotide      Use dinculeotide shuffled z-scores (default=off)");
  printf("%s\n","  -l, --locarnate         Use decision model for structural alignments (default=off)");
  printf("%s\n\n","  -n, --no-shuffle      Do not use shuffling (default=off)");

}

PRIVATE void version(void){
  /*printf("RNAz version " PACKAGE_VERSION "\n");*/
  printf("RNAz version 2.0pre\n");
  exit(EXIT_SUCCESS);
}







