/*********************************************************************
 *                                                                   *
 *                          strand.c                                 *
 *                                                                   *
 *   Functions to predict the reading direction of a                 * 
 *   structured RNA candidate.                                       *
 *                                                                   *
 *                c Kristin Missal                                   *
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
#include "strand.h"

/* wash */
#include "strand.inc"

/* Global variables */
struct svm_model *strand_model=NULL;
double scale_strand_decision[][3]= {    // Boundaries for strand descriptors
{ 0, -0.6, 0.47 },
{ 1, -47.9, 48.35 },
{ 2, -64.44, 66.86 },
{ 3, -5.38, 5.38 },
{ 4, 2, 6 },
{ 5, 35.04, 99.70999999999999 },
{ -1, 0, 0}                            // To check if scale_strand_decision has the same number of entries as svm_node
};


char *strandEnVar = "";        // Name of environment variable which points to the models. Used for easy error printing.

/* PURPOSE: Sets the predictor variables  
   PARAMS:  char *window:            Window of alignment for which strand is predicted
            char *window_reverse:    Reverse complement of window
	    int length:              Length of alignment window
	    int n_seq:               Number of sequences in alignment window 
	    double *deltaSCI:        1st predictor
	    double *deltaMeanMFE:    2nd predictor
	    double *deltaConsMFE:    3rd predictor
	    double *deltaZ:          4th predictor
	    char *modelDir:          Directory with SVM models
   RETURNS: 0 on success
*/
void get_strand_predictors(
    char *window[], 
    char *window_reverse[], 
    int length,
    int n_seq,
    double *deltaSCI, 
    double *deltaMeanMFE,
    double *deltaConsMFE,
    double *deltaZ,
    char *modelDir) 
{
    int i,j;
    char *structure=NULL, *structure_reverse=NULL;
    char *singleStruc, *woGapsSeq, *singleStruc_reverse, *woGapsSeq_reverse;
    double decValue=0;
    double singleMFE,sumMFE,singleZ,sumZ,z,sci;
    double singleMFE_reverse,sumMFE_reverse,singleZ_reverse,sumZ_reverse,z_reverse,sci_reverse;
    double min_en, real_en;
    double min_en_reverse, real_en_reverse;
       
    sumZ=0; sumMFE=0, sumZ_reverse=0, sumMFE_reverse=0;

    /* Init regression SVM for mean and standard deviation of MFE */
    regression_svm_init(modelDir);

    structure = (char *) space((unsigned) length+1);
    structure_reverse = (char *) space((unsigned) length+1);
    min_en = alifold(window, structure);
    min_en_reverse = alifold(window_reverse, structure_reverse);

    /* get z-score, mean mfe, average single mfe and sci*/
    for (i=0;i<n_seq;i++)
    {
	singleStruc = space(strlen(window[i])+1);
	woGapsSeq = space(strlen(window[i])+1);
	singleStruc_reverse = space(strlen(window_reverse[i])+1);
	woGapsSeq_reverse = space(strlen(window_reverse[i])+1);
	j=0;
	while (window[i][j])
	{
	    window[i][j]=toupper(window[i][j]);
	    if (window[i][j]=='T') window[i][j]='U';
	    if (window[i][j]!='-')
	    {
		woGapsSeq[strlen(woGapsSeq)]=window[i][j];
		woGapsSeq[strlen(woGapsSeq)]='\0';
	    }
	    ++j;
	}
	j=0;
	while (window_reverse[i][j])
	{
	    window_reverse[i][j]=toupper(window_reverse[i][j]);
	    if (window_reverse[i][j]=='T') window_reverse[i][j]='U';
	    if (window_reverse[i][j]!='-')
	    {
		woGapsSeq_reverse[strlen(woGapsSeq_reverse)]=window_reverse[i][j];
		woGapsSeq_reverse[strlen(woGapsSeq_reverse)]='\0';
	    }
	    ++j;
	}
	singleMFE = fold(woGapsSeq, singleStruc);
	singleZ=mfe_zscore(woGapsSeq,singleMFE);
	singleMFE_reverse = fold(woGapsSeq_reverse, singleStruc_reverse);
	singleZ_reverse=mfe_zscore(woGapsSeq_reverse,singleMFE_reverse);
	sumZ+=singleZ;
	sumMFE+=singleMFE;
	sumZ_reverse+=singleZ_reverse;
	sumMFE_reverse+=singleMFE_reverse;
	free(woGapsSeq);
	free(singleStruc);
	free(woGapsSeq_reverse);
	free(singleStruc_reverse);
    }

    {
	int i; double s=0;
	extern int eos_debug;
	eos_debug=-1; /* shut off warnings about nonstandard pairs */
	for (i=0; window[i]!=NULL; i++) 
	    s += energy_of_struct(window[i], structure);
	real_en = s/i;
    }
    
    z=sumZ/n_seq;
    sci=min_en/(sumMFE/n_seq);
    z_reverse=sumZ_reverse/n_seq;
    sci_reverse=min_en_reverse/(sumMFE_reverse/n_seq);
    
    *deltaSCI = sci - sci_reverse;
    *deltaMeanMFE = sumMFE/n_seq - sumMFE_reverse/n_seq;
    *deltaConsMFE = min_en - min_en_reverse;
    *deltaZ = z - z_reverse;

    /* Free svm regression model */
    regression_svm_free();

    return;
}

/* PURPOSE: Predicts strand (reading direction) of structured RNA candidate. 
   PARAMS:  double deltaSCI, double deltaMeanMFE, double deltaConsMFE, double deltaZ 
            int n_seq:               Number of sequences in alignment 
	    double id:               Mean pairwise identity of alignment 
            int *strand:             1 if given window is the predicted direction, 0 if reverse complement of given window is predicted 
            double *prob:            Decision probability of SVM that given window is the predicted direction
	    doube *decValue:         Decision value of SVM
	    char *modelDir:          Directory which contains 'strand.model'
   RETURNS: 0 if range error / 1 on success (wash)
*/
int predict_strand(double deltaSCI, double deltaMeanMFE, double deltaConsMFE, double deltaZ, int n_seq, double id, 
	    int *strand, double *prob, double *decValue, char *modelDir)
{
    char outOfRange[1024]="";

	/*wash*/
    /* if (modelDir==NULL) */
/*     { */
/* 	fprintf(stderr, "\n\nERROR: Could not find decision-model 'strand.model' for strand prediction. You have to set the %s" */
/* 		" enviroment variable pointing to the model files!\n\n", strandEnVar); */
/* 	exit(EXIT_FAILURE); */
/*     } */

    strand_svm_init(modelDir);
    
    if (strand_model==NULL)
    {
	fprintf(stderr, "\n\nERROR: Could not find decision-model 'strand.model' for strand prediction. You have to set the %s"
		" enviroment variable pointing to the model files!\n\n", strandEnVar);
	exit(EXIT_FAILURE);
    }
    
    /* check if calculated descriptors are in the models' range*/
    if(checkStrandDescriptors(deltaSCI, deltaMeanMFE, deltaConsMFE, deltaZ, outOfRange) == 1) 
    {
	  //fprintf(stderr, "WARNING: %s", outOfRange);

	  return 0;
	//exit(EXIT_SUCCESS);
    }
    
    /* do classification */
    predict_strand_SCI_MeanMFE_ConsMFE_Z(prob, decValue, deltaSCI, deltaMeanMFE, deltaConsMFE, deltaZ, n_seq, id);
    //predict_strand_SCI_MeanMFE_ConsMFE(prob, decValue, deltaSCI, deltaMeanMFE, deltaConsMFE, n_seq, id);
    //predict_strand_SCI_MeanMFE(prob, decValue, deltaSCI, deltaMeanMFE, n_seq, id);
    //predict_strand_SCI_ConsMFE(prob, decValue, deltaSCI, deltaConsMFE, n_seq, id);
    //predict_strand_SCI_ConsMFE_Z(prob, decValue, deltaSCI, deltaConsMFE, deltaZ, n_seq, id);
    //predict_strand_SCI_Z(prob, decValue, deltaSCI, deltaZ, n_seq, id);
    //predict_strand_ConsMFE_Z(prob, decValue, deltaConsMFE, deltaZ, n_seq, id);
    //predict_strand_MEANMFE_CONSMFE(prob, decValue, deltaMeanMFE, deltaConsMFE, n_seq, id);
    //predict_strand_MEANMFE_CONSMFE_Z(prob, decValue, deltaMeanMFE, deltaConsMFE, deltaZ, n_seq, id);

    /* Check if probability is a true probability */
    if(*prob<0 || *prob>1 )
    {
	fprintf(stderr, "\n\nERROR: SVM classification probability is %f!\n\n", *prob);
	exit(EXIT_SUCCESS);
    }
  
    if(*prob > 0.5)
	*strand = 1;
    else
	*strand = 0;

    strand_svm_free();

    return 1;
}

/*
  PURPOSE: Check if descriptors are within range
  PARAMS:  double deltaSCI, double deltaMeanMFE, double deltaConsMFE, char* error
  RETURNS: 1 if descriptors are within range, else 0
*/
int checkStrandDescriptors(double deltaSCI, double deltaMeanMFE, double deltaConsMFE, double deltaZ, char *error) {
	
  if(deltaSCI > scale_strand_decision[0][2] || deltaSCI < scale_strand_decision[0][1]) 
    {
	strcpy(error, " deltaSCI is out of range ");	
    }					
  if(deltaMeanMFE > scale_strand_decision[1][2] || deltaMeanMFE < scale_strand_decision[1][1]) 
    {
	strcpy(error, " deltaMeanMFE is out of range ");	
    }					
  if(deltaConsMFE > scale_strand_decision[2][2] || deltaConsMFE < scale_strand_decision[2][1]) 
    {
	strcpy(error, " deltaConsMFE is out of range\n");	
    }
  if(deltaZ > scale_strand_decision[3][2] || deltaZ < scale_strand_decision[3][1]) 
    {
	strcpy(error, " deltaZ is out of range\n");	
    }
  
  if(strcmp(error, "") != 0) return 1;
  else return 0;
}


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

/* wash */
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
void get_strand_model(char *basefilename){

    char fn[256];
    
    if (basefilename!=NULL)
    {
	strncpy(fn, basefilename, 238); strcat(fn, "strand.model");
	strand_model = svm_load_model(fn);
	//if (strand_model!=NULL) 
	//    fprintf(stderr, "Loaded: %s (%p)\n",fn,strand_model);    
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
void strand_svm_free()
{
    svm_destroy_model(strand_model);
    strand_model=NULL;
    return;
}


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

void strand_env_variable_init(char* envVariable)
{
    strandEnVar = envVariable;
}

/****************************************************************************************************************
  Several functions which do prediction based on different combinations of classification variables
*****************************************************************************************************************/
void predict_strand_SCI_MeanMFE_ConsMFE(double* prob, double* decValue, double deltaSCI, double deltaMeanMFE, double deltaConsMFE, int n_seq, double id)
{
  struct svm_node node[6];
  double value=0;
 
  node[0].index = 1; node[0].value = deltaSCI;
  node[1].index = 2; node[1].value = deltaMeanMFE;
  node[2].index = 3; node[2].value = deltaConsMFE;
  node[3].index = 4; node[3].value = (double)n_seq;
  node[4].index = 5; node[4].value = id;
  node[5].index = -1;

  scale_strand_decision_node((struct svm_node*)&node);

  svm_predict_values(strand_model,node,&value);
  *decValue=value;
  svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                      // to be C_SVC or nu_SVC!
  *prob=value;
 
  return;
}

void predict_strand_SCI_MeanMFE_ConsMFE_Z(double* prob, double* decValue, double deltaSCI, double deltaMeanMFE, 
					  double deltaConsMFE, double deltaZ, int n_seq, double id)
{
  struct svm_node node[7];
  double value=0;
 
  node[0].index = 1; node[0].value = deltaSCI;
  node[1].index = 2; node[1].value = deltaMeanMFE;
  node[2].index = 3; node[2].value = deltaConsMFE;
  node[3].index = 4; node[3].value = deltaZ;
  node[4].index = 5; node[4].value = (double)n_seq;
  node[5].index = 6; node[5].value = id;
  node[6].index = -1;

  scale_strand_decision_node((struct svm_node*)&node);

  svm_predict_values(strand_model,node,&value);
  *decValue=value;
  svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                      // to be C_SVC or nu_SVC!
  *prob=value;
 
  return;
}

void predict_strand_SCI_MeanMFE(double* prob, double* decValue, double deltaSCI, double deltaMeanMFE, int n_seq, double id)
{
  struct svm_node node[5];
  double value=0;
 
  node[0].index = 1; node[0].value = deltaSCI;
  node[1].index = 2; node[1].value = deltaMeanMFE;
  node[2].index = 3; node[2].value = (double)n_seq;
  node[3].index = 4; node[3].value = id;
  node[4].index = -1;

  scale_strand_decision_node((struct svm_node*)&node);
  
  svm_predict_values(strand_model,node,&value);
  *decValue=value;
  svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                      // to be C_SVC or nu_SVC!
  *prob=value;
 
  return;
}

void predict_strand_SCI_ConsMFE(double* prob, double* decValue, double deltaSCI, double deltaConsMFE, int n_seq, double id)
{
  struct svm_node node[5];
  double value=0;
 
  node[0].index = 1; node[0].value = deltaSCI;
  node[1].index = 2; node[1].value = deltaConsMFE;
  node[2].index = 3; node[2].value = (double)n_seq;
  node[3].index = 4; node[3].value = id;
  node[4].index = -1;

  scale_strand_decision_node((struct svm_node*)&node);

  svm_predict_values(strand_model,node,&value);
  *decValue=value;
  svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                      // to be C_SVC or nu_SVC!
  *prob=value;
 
  return;
}

void predict_strand_SCI_Z(double* prob, double* decValue, double deltaSCI, double deltaZ, int n_seq, double id)
{
  struct svm_node node[5];
  double value=0;
 
  node[0].index = 1; node[0].value = deltaSCI;
  node[1].index = 2; node[1].value = deltaZ;
  node[2].index = 3; node[2].value = (double)n_seq;
  node[3].index = 4; node[3].value = id;
  node[4].index = -1;

  scale_strand_decision_node((struct svm_node*)&node);

  svm_predict_values(strand_model,node,&value);
  *decValue=value;
  svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                      // to be C_SVC or nu_SVC!
  *prob=value;
 
  return;
}

void predict_strand_ConsMFE_Z(double* prob, double* decValue, double deltaConsMFE, double deltaZ, int n_seq, double id)
{
  struct svm_node node[5];
  double value=0;
 
  node[0].index = 1; node[0].value = deltaConsMFE;
  node[1].index = 2; node[1].value = deltaZ;
  node[2].index = 3; node[2].value = (double)n_seq;
  node[3].index = 4; node[3].value = id;
  node[4].index = -1;

  scale_strand_decision_node((struct svm_node*)&node);

  svm_predict_values(strand_model,node,&value);
  *decValue=value;
  svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                      // to be C_SVC or nu_SVC!
  *prob=value;
 
  return;
}

void predict_strand_SCI_ConsMFE_Z(double *prob, double *decValue, double deltaSCI, double deltaConsMFE, double deltaZ, int n_seq, double id)
{
    struct svm_node node[6];
    double value=0;
    
    node[0].index = 1; node[0].value = deltaSCI;
    node[1].index = 2; node[1].value = deltaConsMFE;
    node[2].index = 3; node[2].value = deltaZ;
    node[3].index = 4; node[3].value = (double)n_seq;
    node[4].index = 5; node[4].value = id;
    node[5].index = -1;
    
    scale_strand_decision_node((struct svm_node*)&node);
    
    svm_predict_values(strand_model,node,&value);
    *decValue=value;
    svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                        // to be C_SVC or nu_SVC!
    *prob=value;
    
    return;

}


void predict_strand_MEANMFE_CONSMFE(double *prob, double *decValue, double deltaMeanMFE, double deltaConsMFE, int n_seq, double id)
{
    struct svm_node node[5];
    double value=0;
    
    node[0].index = 1; node[0].value = deltaMeanMFE;
    node[1].index = 2; node[1].value = deltaConsMFE;
    node[2].index = 3; node[2].value = (double)n_seq;
    node[3].index = 4; node[3].value = id;
    node[4].index = -1;
    
    scale_strand_decision_node((struct svm_node*)&node);
    
    svm_predict_values(strand_model,node,&value);
    *decValue=value;
    svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                        // to be C_SVC or nu_SVC!
    *prob=value;
    
    return;
}

void predict_strand_MEANMFE_CONSMFE_Z(double *prob, double *decValue, double deltaMeanMFE, double deltaConsMFE, double deltaZ, int n_seq, double id)
{
    struct svm_node node[6];
    double value=0;
    
    node[0].index = 1; node[0].value = deltaMeanMFE;
    node[1].index = 2; node[1].value = deltaConsMFE;
    node[2].index = 3; node[2].value = deltaZ;
    node[3].index = 4; node[3].value = (double)n_seq;
    node[4].index = 5; node[4].value = id;
    node[5].index = -1;
    
    scale_strand_decision_node((struct svm_node*)&node);
    
    svm_predict_values(strand_model,node,&value);
    *decValue=value;
    svm_predict_probability(strand_model,node,&value);  // ATTENTION: If model parameters are not set correctly, then value is not set, type of SVM need to 
                                                        // to be C_SVC or nu_SVC!
    *prob=value;
    
    return;
}




