/*********************************************************************
 *                                                                   *
 *                          strand.c                                 *
 *                                                                   *
 *   Functions to predict the reading direction of a                 * 
 *   structured RNA candidate.                                       *
 *                                                                   *
 *                c Kristin Reiche                                   *
 *   $Id: strand.c,v 1.2 2008-01-24 10:20:10 wash Exp $                                                          *
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
#include "strand.h"

#define MAX_DESCRIPTORS 8

/* wash */
#include "strand.inc"

/* Global variables */
static struct svm_model *strand_model=NULL;

/*wash*/

/*static double scale_strand_decision[MAX_DESCRIPTORS+1][3];*/

double scale_strand_decision[][3]= {    // Boundaries for strand descriptors
  { 0, -0.5600000000000001,0.48 },
  { 1,-4.92, 4.92 },
  { 2, -48.01, 47.64},
  { 3, -66.78, 64.44},
  { 4,  2, 6},
  { 5,  46.74, 99.54000000000001},
  { 6, 0, 48.53},
  { -1, 0, 0}                            // To check if scale_strand_decision has the same number of entries as svm_node
};



static char *strandEnVar = "";        		// Name of environment variable which points to the models. Used for easy error printing.
static double fSCI, fMeanMFE, fConsMFE, fZ, rSCI, rMeanMFE, rConsMFE, rZ;  // global variables which save the values for each strand
static char *output;				// Holds sequences, structures and their MFEs for printing
static char outOfRangeWarning[1024] = "";            // If descriptors are not within training range, error is reported in outOfRange

/* PURPOSE: Sets the predictor variables  
   PARAMS:  char *window:            Window of alignment for which strand is predicted
            char *window_reverse:    Reverse complement of window
	    char *names:	     Sequences names
	    int length:              Length of alignment window
            int length_rev:          Length of reverse alignment window
	    int n_seq:               Number of sequences in alignment window 
	    double *deltaSCI:        1st predictor
	    double *deltaMeanMFE:    2nd predictor
	    double *deltaConsMFE:    3rd predictor
	    double *deltaZ:          4th predictor
	    double *freqGU:          Prediction depends on GU base pair content in alignment.
	    char *modelDir:          Directory with SVM models
   RETURNS: -
*/
void get_strand_predictors(
    char *window[], 
    char *window_reverse[], 
    char *names[],
    int length,
    int length_rev,
    int n_seq,
    double *deltaSCI, 
    double *deltaMeanMFE,
    double *deltaConsMFE,
    double *deltaZ,
    double *consFreqGU,
    char *modelDir) 
{
    int i,j;
    char *structure=NULL, *structure_reverse=NULL;
    char *singleStruc=NULL, *woGapsSeq=NULL, *singleStruc_reverse=NULL, *woGapsSeq_reverse=NULL;
    char *string=NULL, *o=NULL, *o_reverse=NULL;
    double decValue=0;
    double singleMFE,sumMFE,singleZ,sumZ,z,sci;
    double singleMFE_reverse,sumMFE_reverse,singleZ_reverse,sumZ_reverse,z_reverse,sci_reverse;
    double min_en, real_en;
    double min_en_reverse, real_en_reverse;
    double sumFreqGU=0.0;   
 
    /* Global RNA package variables */
    dangles=2;
    sumZ=0; sumMFE=0, sumZ_reverse=0, sumMFE_reverse=0;

    output=(char *)space(sizeof(char)*(length+16)*(n_seq+1)*3*2);  // Multiply with 2 as positive and negative strands are regarded
    o=(char *)space(sizeof(char)*(length+16)*(n_seq+1)*3);
    o_reverse=(char *)space(sizeof(char)*(length_rev+16)*(n_seq+1)*3);

    /* Replace any T's to U's before calling alifold and fold */
    for (i=0;i<n_seq;i++){
      j=0;
      while (window[i][j]){
        window[i][j]=toupper(window[i][j]);
        if (window[i][j]=='T') window[i][j]='U';
        ++j;
      }
      j=0;
      while (window_reverse[i][j]){
        window_reverse[i][j]=toupper(window_reverse[i][j]);
        if (window_reverse[i][j]=='T') window_reverse[i][j]='U';
        ++j;
      }
    }
    
    /* Init regression SVM for mean and standard deviation of MFE */
    regression_svm_init(modelDir);
    structure = (char *) space((unsigned) length+1);
    structure_reverse = (char *) space((unsigned) length_rev+1);
    min_en = alifold(window, structure);
    free_alifold_arrays();  
    min_en_reverse = alifold(window_reverse, structure_reverse);
    free_alifold_arrays();

    /* Get z-score, mean mfe, average single mfe and sci for forward and reverse window. 
       Get frequency of GU pairs according to number of all base pairs for forward window.
    */
    for (i=0;i<n_seq;i++){
      singleStruc = space(strlen(window[i])+1);
      woGapsSeq = space(strlen(window[i])+1);
      singleStruc_reverse = space(strlen(window_reverse[i])+1);
      woGapsSeq_reverse = space(strlen(window_reverse[i])+1);
      j=0;
      while (window[i][j]){
        if (window[i][j]!='-'){
          woGapsSeq[strlen(woGapsSeq)]=window[i][j];
          woGapsSeq[strlen(woGapsSeq)]='\0';
        }
        ++j;
      }
      j=0;
      while (window_reverse[i][j]){
        if (window_reverse[i][j]!='-'){
          woGapsSeq_reverse[strlen(woGapsSeq_reverse)]=window_reverse[i][j];
          woGapsSeq_reverse[strlen(woGapsSeq_reverse)]='\0';
        }
        ++j;
      }
      singleMFE = fold(woGapsSeq, singleStruc);
      free_arrays();   
      singleZ=mfe_zscore(woGapsSeq,singleMFE);
      singleMFE_reverse = fold(woGapsSeq_reverse, singleStruc_reverse);
      free_arrays();
      singleZ_reverse = mfe_zscore(woGapsSeq_reverse,singleMFE_reverse);
      
      sumZ+=singleZ;
      sumMFE+=singleMFE;
      sumZ_reverse+=singleZ_reverse;
      sumMFE_reverse+=singleMFE_reverse;
      
      sprintf(o+strlen(o),">%s +\n%s\n%s ( %6.2f)\n",
              names[i],woGapsSeq,singleStruc,singleMFE);
      sprintf(o_reverse+strlen(o_reverse),">%s -\n%s\n%s ( %6.2f)\n",
              names[i],woGapsSeq_reverse,singleStruc_reverse,singleMFE_reverse);
      
      free(woGapsSeq);
      free(singleStruc);
      free(woGapsSeq_reverse);
      free(singleStruc_reverse);
    }
    /* consensus structure */
    {
      int i; double s=0;
      extern int eos_debug;
      eos_debug=-1; /* shut off warnings about nonstandard pairs */
      for (i=0; window[i]!=NULL; i++) 
        s += energy_of_struct(window[i], structure);
      real_en = s/i;
      s=0;
      for (i=0; window_reverse[i]!=NULL; i++)
        s += energy_of_struct(window_reverse[i], structure_reverse);
      real_en_reverse = s/i;
    }
    string = (char*)strand_consensus((const char **) window);
    sprintf(o+strlen(o),
            ">consensus +\n%s\n%s (%6.2f = %6.2f + %6.2f) \n",
            string, structure, min_en, real_en, min_en-real_en );   
    
    sumFreqGU = (double)frequency_gu_pairs(string, structure);
    
    string = (char*)strand_consensus((const char **) window_reverse);
    sprintf(o_reverse+strlen(o_reverse),
            ">consensus -\n%s\n%s (%6.2f = %6.2f + %6.2f) \n",
            string, structure_reverse, min_en_reverse, real_en_reverse, min_en_reverse-real_en_reverse );
    sprintf(output+strlen(output),"%s\n%s",o, o_reverse);
    
    sumFreqGU += (double)frequency_gu_pairs(string, structure_reverse);
    
    z=sumZ/n_seq;
    sci=min_en/(sumMFE/n_seq);
    z_reverse=sumZ_reverse/n_seq;
    sci_reverse=min_en_reverse/(sumMFE_reverse/n_seq);
    
    *deltaSCI = sci - sci_reverse;
    *deltaMeanMFE = sumMFE/n_seq - sumMFE_reverse/n_seq;
    *deltaConsMFE = min_en - min_en_reverse;
    *deltaZ = z - z_reverse;
    *consFreqGU = (double)sumFreqGU; // /2.0;

    // Save values in global variables
    fSCI = sci;
    fMeanMFE = sumMFE/n_seq;
    fConsMFE = min_en;
    fZ = z;
    rSCI = sci_reverse;
    rMeanMFE = sumMFE_reverse/n_seq;
    rConsMFE = min_en_reverse;
    rZ = z_reverse;

    /* Free svm regression model */
    regression_svm_free();

    free(o); free(o_reverse);
    free(structure); free(structure_reverse);
    free(string);
    return;
}

/* PURPOSE: Predicts strand (reading direction) of structured RNA candidate. 
   PARAMS:  double deltaSCI, double deltaMeanMFE, double deltaConsMFE, double deltaZ 
            int n_seq:               Number of sequences in alignment 
	    double id:               Mean pairwise identity of alignment: (mpi(+)+mpi(-))/2
            double consFreqGU:       GU base pair frequency in consensus structures
	    double *score:           RNAstrand score
            double *prob:            Decision probability of SVM for positive reading direction
	    doube *decValue:         Decision value of SVM
	    char *modelDir:          Directory which contains 'strand.model'
            int code: 		     Code which combination of classification variables should be taken
				     (usually: CODE_SCI_Z_MEANMFE_CONSMFE)
   RETURNS: -
*/
int strand(double deltaSCI, double deltaMeanMFE, double deltaConsMFE, double deltaZ, int n_seq, double id, double consFreqGU, 
	    double *score, double *prob, double *decValue, char *modelDir, int code)
{
    double descriptors[MAX_DESCRIPTORS];
    FILE* scale_file;
    char fn[256];
    int n_descr=0;	// Number of chosen descriptors
    
    /*wash*/
    /*if (modelDir==NULL)
      {
        fprintf(stderr, "\n\nERROR: Could not find decision-model 'strand.model' for strand prediction. You have to set the %s"
                " enviroment variable pointing to the model files!\n\n", strandEnVar);
        exit(EXIT_FAILURE);
        }*/
    
    strand_svm_init(modelDir);

    /* Get file containing descriptor ranges. */
    if (modelDir!=NULL){
        strncpy(fn, modelDir, 238); strcat(fn, "/strand.range");
        scale_file = fopen(fn, "r");
                if(scale_file == NULL) {
                  	fprintf(stderr, "\n\nERROR: Could not find descriptor ranges 'strand.range'. You have to set the %s"
                		" enviroment variable pointing to the model files!\n\n", strandEnVar);
			exit(EXIT_FAILURE);
		}
    }

    /*wash*/
    /* Ranges to global vector scale_strand_decision */
    /*read_scale_vector(scale_file); */  

    
    if (strand_model==NULL)
      {
        fprintf(stderr, "\n\nERROR: Could not find decision-model 'strand.model' for strand prediction. You have to set the %s"
                " enviroment variable pointing to the model files!\n\n", strandEnVar);
        exit(EXIT_FAILURE);
      }
    
    /* get classification variables */
    if(code==CODE_SCI){
    	descriptors[0] = deltaSCI;
	n_descr=1;
    }
    else if(code==CODE_MEANMFE){ 
	descriptors[0] = deltaMeanMFE;
    	n_descr=1;
    }
    else if(code==CODE_CONSMFE){
    	descriptors[0] = deltaConsMFE;
	n_descr=1;
    }
    else if(code==CODE_Z){
    	descriptors[0] = deltaZ;
	n_descr=1;
    }
    else if(code==CODE_SCI_Z){
	descriptors[0] = deltaSCI;
	descriptors[1] = deltaZ;
	n_descr=2;
    }
    else if(code==CODE_SCI_Z_MEANMFE){
	descriptors[0] = deltaSCI;
	descriptors[1] = deltaZ;
	descriptors[2] = deltaMeanMFE;
	n_descr=3;
    }
    else if(code==CODE_SCI_Z_CONSMFE){
	descriptors[0] = deltaSCI;
	descriptors[1] = deltaZ;
	descriptors[2] = deltaConsMFE;
	n_descr=3;
    }
    else if(code==CODE_SCI_Z_MEANMFE_CONSMFE){
	descriptors[0] = deltaSCI;
	descriptors[1] = deltaZ;
	descriptors[2] = deltaMeanMFE;
	descriptors[3] = deltaConsMFE;
	n_descr=4;
    }
    else if(code==CODE_MEANMFE_CONSMFE){
	descriptors[0] = deltaMeanMFE;
        descriptors[1] = deltaConsMFE;
        n_descr=2;
    } 
    else if(code==CODE_Z_MEANMFE_CONSMFE){
	descriptors[0] = deltaZ;
        descriptors[1] = deltaMeanMFE;
        descriptors[2] = deltaConsMFE;
        n_descr=3;
    }
    else if(code==CODE_Z_CONSMFE){
	descriptors[0] = deltaZ;
        descriptors[1] = deltaConsMFE;
        n_descr=2;
    }
    else if(code == CODE_SCI_CONSMFE){
      descriptors[0] = deltaSCI;
      descriptors[1] = deltaConsMFE;
      n_descr=2;
    }
    else if(code == CODE_SCI_MEANMFE){
      descriptors[0] = deltaSCI;
      descriptors[1] = deltaMeanMFE;
      n_descr=2;
    }
    else if(code == CODE_SCI_MEANMFE_CONSMFE){
      descriptors[0] = deltaSCI;
      descriptors[1] = deltaMeanMFE;
      descriptors[2] = deltaConsMFE;
      n_descr=3;
    }
    else if(code == CODE_Z_MEANMFE){
      descriptors[0] = deltaZ;
      descriptors[1] = deltaMeanMFE;
      n_descr=2;
    }

    /* check if descriptors are in the models' range. If not, warning is printed to outOfRangeWarning. */
    if(checkStrandDescriptors(descriptors, n_descr, outOfRangeWarning, code) == 1){
        //do nothing: Error message is printed in checkStrandDescriptors
      return 0;
    }
    /* do classification */
    predict_strand(prob, decValue, descriptors, n_seq, id, consFreqGU, n_descr);

    /* Check if probability is a true probability */
    if(*prob<0 || *prob>1 )
    {
	fprintf(stderr, "\n\nERROR: SVM classification probability is %f!\n\n", *prob);
	exit(EXIT_SUCCESS);
    }
  
    *score = 2. * (*prob) - 1.;

    strand_svm_free();

    return 1;
}


/*  
   PURPOSE: Calculates frequency of GU base pairs of a given structure.
   PARAMS:  char *string: sequence
            char *struc:  structure
   RETURNS: Frequency 
*/
double frequency_gu_pairs(char *string, char *struc)
{
  int i, length, cntGU=0, cntAll=0, pos=0;
  char *stack=NULL;
  char *pair=(char *)space(sizeof(char)*2);

  length = (int) strlen(string);
  if(length != (int) strlen(struc)){
        fprintf(stderr, "\n\nERROR: Sequence and structure have different length (seq: %d struc: %d)! \n\n", strlen(string), strlen(struc));
        exit(EXIT_FAILURE);
  }

  stack=(char *)space(sizeof(char)*((length/2)+1));

  for(i=0; i<length; i++){
        if(struc[i]=='('){
                stack[pos]=string[i];
                pos++;
        }
        if(struc[i]==')'){
                pair[0] = toupper(stack[pos-1]);
                pair[1] = toupper(string[i]);
                pos--;
                if((pair[0]=='G' && pair[1]=='U')
                        || (pair[0]=='U' && pair[1]=='G')
                        || (pair[0]=='G' && pair[1]=='T')
                        || (pair[0]=='T' && pair[1]=='G')
                        ){
                        cntGU++; cntAll++;
                }
                else{ cntAll++; }
        }
  }
  if(cntAll==0){
        return 0;
  }
  return((double)((double)cntGU/(double)cntAll)*100.0);
  
}

/* Calculates consensus of alignment 
*  Code taken from RNAz 0.1.1 */
static char *strand_consensus(const char *AS[]) 
{
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


/****************************************************************************************************************
  Read and scale descriptor values.
*****************************************************************************************************************/
/*
   PURPOSE: Read unscaled range of descriptors from file strand.range to global vector scale_strand_decision. 
   PARAMS:  Input file 'strand.range'.
   RETURNS: Number of descriptors
*/
static int read_scale_vector(FILE *scaleFile)
{

  char *line;
  int i=0, n_descriptors=0;

  if ((line=(char *)get_line(scaleFile)) == NULL) {
        fprintf(stderr, "ERROR: Empty scale file\n"); return 0;
  }
  while (line!=NULL) {
        scale_strand_decision[i][0]=i;
        sscanf(line,"%lf %lf", &scale_strand_decision[i][1], &scale_strand_decision[i][2]);
        free(line);
        //printf("'%f %f'\n", scale_strand_decision[i][1], scale_strand_decision[i][2]);
        line = (char *)get_line(scaleFile);
        ++i;
        ++n_descriptors;
  }
  /* To check if scale_strand_decision has the same number of entries as svm_node. */
  scale_strand_decision[i][0]=-1; scale_strand_decision[i][1]=0; scale_strand_decision[i][2]=0;

  return n_descriptors;
}

/*
  PURPOSE: Check if descriptors are within range
  PARAMS:  double deltaSCI, double deltaMeanMFE, double deltaConsMFE, char* error
           int code: Codes used combination of descriptors. This gives us the descriptor which is out of range.
  RETURNS: 1 if descriptors are within range, else 0
*/
int checkStrandDescriptors(double* descriptors, int n_descr, char *error, int code) 
{
  int i;
  char descr[1024]="";	
  for(i=0; i<n_descr; i++){
	if(descriptors[i] > scale_strand_decision[i][2] || descriptors[i] < scale_strand_decision[i][1]){
		i2descriptor(code, i, descr);
		sprintf(error, "Descriptor '%s' is out of range (min=%.4lf max=%.4lf).", descr, scale_strand_decision[i][1], scale_strand_decision[i][2]);
    	}
  }
  if(strcmp(error, "") != 0) return 1;
  else return 0;
}

/* original code taken from RNAz 0.1.1 */
void scale_strand_node_linear(struct svm_node* node,double scale[][3])
{
    int i=0;
    double min,max;
    double from,to;
    from=-1;
    to=1;

    // Correct number of nodes as given in array scale?
    int cntNode=0;
    while(1)
    {
        cntNode++;
        ++i;
        if (node[i].index==-1) break;
    }
    i=0;
    int cntScale=0;
    while(1)
    {
        cntScale++;
        ++i;
        if(scale[i][0]==-1) break;
    }
    i=0;
    if(cntNode != cntScale)
    {
        fprintf(stderr, "cntNode=%d cntScale=%d", cntNode, cntScale);
        nrerror("ERROR: Number of svm nodes does not correspond to number of entries found in scale array (scale_strand_node_linear)\n\n");
    }

    while (1)
    {
        min=scale[node[i].index-1][1];
        max=scale[node[i].index-1][2];
        node[i].value=from+(to-from)*(node[i].value-min)/(max-min);
        ++i;
        if (node[i].index==-1) break;
    }
}


void scale_strand_decision_node(struct svm_node* node)
{
    scale_strand_node_linear(node,scale_strand_decision);
}


/****************************************************************************************************************
  Init, read and free SVM model.
*****************************************************************************************************************/
/*
  PURPOSE: Initialization of svm model, which predicts strand of structured ncRNA.
  PARAMS:  char* basefilename: Name of model directory
  RETURNS: -
*/

/* void strand_svm_init(char *basefilename) */
/* { */
/*     strand_model=NULL; */
    
/*     get_strand_model(basefilename); */
/*     if (strand_model==NULL) */
/*     { */
/*         fprintf(stderr, "\n\nERROR: Could not load strand model 'strand.model'. \n" */
/* 		"You have to set the %s enviroment variable " */
/* 		"pointing to the model files!\n\n", strandEnVar); */
/* 	exit(EXIT_FAILURE); */
/*     } */
/*     return; */
/* } */

/*wash*/
void strand_svm_init(char *basefilename){
  strand_model=NULL;
  
  strand_model=svm_load_model_string(default_strand_model_string);
  
  if (strand_model==NULL){
    fprintf(stderr, "\n\nERROR: Could not load strand model. \n");
    exit(EXIT_FAILURE);
  }
  return;
}








/*
  PURPOSE: Sets the svm model
  PARAMS:  char* basefilename: Name of model directory
  RETURNS: -
*/
void get_strand_model(char *basefilename)
{

    char fn[256];

    if (basefilename!=NULL)
    {
	strncpy(fn, basefilename, 238); strcat(fn, "/strand.model");
	strand_model = svm_load_model(fn);
	if (strand_model==NULL){ 
	    fprintf(stderr, "\n\nERROR:Could not find decision-model for strand prediction. You have to set "
                "the %s enviroment variable pointing to the model files!\n\n", strandEnVar);    
	    exit(EXIT_FAILURE);
	}
    } 
    else 
    {
	fprintf(stderr, "\n\nERROR: Could not find decision-model for strand prediction. You have to set "
		"the %s enviroment variable pointing to the model files!\n\n", strandEnVar);
	exit(EXIT_FAILURE);
    }

    return;
}


/* 
   PURPOSE: Destroys SVM-models 
   PARAMS:  -
   RETURNS: -
*/
static void strand_svm_free()
{
    svm_destroy_model(strand_model);
    strand_model=NULL;
    return;
}


/****************************************************************************************************************
 Prediction  
*****************************************************************************************************************/
void predict_strand(double* prob, double* decValue, double* descriptors, int n_seq, double id, double consFreqGU, int n_descr)
{
  struct svm_node node[MAX_DESCRIPTORS];
  double value[2]={0,0};
  int i=0;

  for(i=0; i<n_descr; i++){
  	node[i].index = i+1; node[i].value = descriptors[i];
  } 
  node[i].index = i+1; node[i].value = (double)n_seq; i++;
  node[i].index = i+1; node[i].value = id; i++;
  node[i].index = i+1; node[i].value = consFreqGU; i++;
  node[i].index = -1;
  scale_strand_decision_node((struct svm_node*)&node);

  svm_predict_values(strand_model,node,value);
  *decValue=value[0];
  svm_predict_probability(strand_model,node,value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                      // to be C_SVC or nu_SVC!
  *prob=value[0];
 
  return;
}

/****************************************************************************************************************
  Functions to get prediction values for each strand
*****************************************************************************************************************/
double get_strand_forward_SCI()
{ return fSCI; }
double get_strand_forward_meanMFE()
{ return fMeanMFE; }
double get_strand_forward_consMFE()
{ return fConsMFE; }
double get_strand_forward_Z()
{ return fZ; }
double get_strand_reverse_SCI()
{ return rSCI; }
double get_strand_reverse_meanMFE()
{ return rMeanMFE; }
double get_strand_reverse_consMFE()
{ return rConsMFE; }
double get_strand_reverse_Z()
{ return rZ; }
char* get_strand_output()
{ return output; }
char* get_strand_warning()
{ return outOfRangeWarning; }

/****************************************************************************************************************
  Functions which free dynamic global variables.  
*****************************************************************************************************************/
void free_strand()
{ free(output); }



/****************************************************************************************************************
  Miscellaneous
*****************************************************************************************************************/
void strand_env_variable_init(char* envVariable)
{
    strandEnVar = envVariable;
}

void code2string(int code)
{
  if(code == CODE_SCI)
    printf("SCI (internal:%d) ", CODE_SCI);
  else if(code == CODE_MEANMFE)
    printf("MeanMFE (internal:%d) ", CODE_MEANMFE);
  else if(code == CODE_CONSMFE)
    printf("ConsMFE (internal:%d) ", CODE_CONSMFE);
  else if(code == CODE_Z)
    printf("Z (internal:%d) ", CODE_Z);
  else if(code == CODE_SCI_Z)
    printf("SCI_Z (internal:%d) ", CODE_SCI_Z);
  else if(code == CODE_SCI_Z_MEANMFE)
    printf("SCI_Z_MeanMFE (internal:%d) ", CODE_SCI_Z_MEANMFE);
  else if(code == CODE_SCI_Z_CONSMFE)
    printf("SCI_Z_ConsMFE (internal:%d) ", CODE_SCI_Z_CONSMFE);
  else if(code == CODE_SCI_Z_MEANMFE_CONSMFE)
    printf("SCI_Z_MeanMFE_ConsMFE (internal:%d) ", CODE_SCI_Z_MEANMFE_CONSMFE);
  else if(code == CODE_MEANMFE_CONSMFE)
    printf("MeanMFE_ConsMFE (internal:%d) ", CODE_MEANMFE_CONSMFE);
  else if(code == CODE_Z_MEANMFE_CONSMFE)
    printf("Z_MeanMFE_ConsMFE (internal:%d) ", CODE_Z_MEANMFE_CONSMFE);
  else if(code == CODE_Z_CONSMFE)
    printf("Z_ConsMFE (internal:%d) ", CODE_Z_CONSMFE);
  else if(code == CODE_SCI_CONSMFE)
    printf("SCI_ConsMFE (internal:%d) ", CODE_SCI_CONSMFE);
  else if(code == CODE_SCI_MEANMFE)
    printf("SCI_MeanMFE (internal:%d) ", CODE_SCI_MEANMFE);
  else if(code == CODE_SCI_MEANMFE_CONSMFE)
    printf("SCI_MeanMFE_ConsMFE (internal:%d) ", CODE_SCI_MEANMFE_CONSMFE);
  else if(code == CODE_Z_MEANMFE)
    printf("Z_MeanMFE (internal:%d)", CODE_Z_MEANMFE);

  return;
}

void i2descriptor(int code, int i, char* descr)
{
   char dsci[1024] = "delta structure conservation index";
   char dmeanmfe[1024] = "delta mean single sequence MFE";
   char dconsmfe[1024] = "delta consensus MFE";
   char dz[1024] = "delta mean z-score";
 
   switch (code){
	case CODE_SCI:
		strcpy(descr, dsci);
        	break;
        case CODE_MEANMFE:
        	strcpy(descr, dmeanmfe);
		break;
        case CODE_CONSMFE:
        	strcpy(descr, dconsmfe);
		break;
	case CODE_Z:
        	strcpy(descr, dz);
		break;
    	case CODE_SCI_Z:
		if(i==0){ strcpy(descr, dsci); };
		if(i==1){ strcpy(descr, dz); }
        	break;
    	case CODE_SCI_Z_MEANMFE:
		if(i==0){ strcpy(descr, dsci); }
		if(i==1){ strcpy(descr, dz); }
		if(i==2){ strcpy(descr, dmeanmfe); }
		break;
    	case CODE_SCI_Z_CONSMFE:
		if(i==0){ strcpy(descr, dsci); }
		if(i==1){ strcpy(descr, dz); }
		if(i==2){ strcpy(descr, dconsmfe); }
		break;
        case CODE_SCI_Z_MEANMFE_CONSMFE:
		if(i==0){ strcpy(descr, dsci); }
		if(i==1){ strcpy(descr, dz); }
		if(i==2){ strcpy(descr, dmeanmfe); }
		if(i==3){ strcpy(descr, dconsmfe); }
		break;
    	case CODE_MEANMFE_CONSMFE:
		if(i==0){ strcpy(descr, dmeanmfe); }
		if(i==1){ strcpy(descr, dconsmfe); }
		break;
	case CODE_Z_MEANMFE_CONSMFE:
		if(i==0){ strcpy(descr, dz); }
		if(i==1){ strcpy(descr, dmeanmfe); }
		if(i==2){ strcpy(descr, dconsmfe); }
		break;
    	case CODE_Z_CONSMFE:
		if(i==0){ strcpy(descr, dz); }
		if(i==1){ strcpy(descr, dconsmfe); }
		break;
    	case CODE_SCI_CONSMFE:
		if(i==0){ strcpy(descr, dsci); }
		if(i==1){ strcpy(descr, dconsmfe); }
		break;
    	case CODE_SCI_MEANMFE:
		if(i==0){ strcpy(descr, dsci); }
		if(i==1){ strcpy(descr, dmeanmfe); }
		break;
    	case CODE_SCI_MEANMFE_CONSMFE:
		if(i==0){ strcpy(descr, dsci); }
		if(i==1){ strcpy(descr, dmeanmfe); }
		if(i==2){ strcpy(descr, dconsmfe); }
		break;
    	case CODE_Z_MEANMFE:
		if(i==0){ strcpy(descr, dz); }
		if(i==1){ strcpy(descr, dmeanmfe); }
		break;
    }
}
