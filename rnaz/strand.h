/*********************************************************************
 *                                                                   *
 *                          strand.h                                 *
 *                                                                   *
 *   Functions to predict the reading direction of a                 *
 *   structured RNA candidate.                                       *
 *                                                                   *
 *                c Kristin Reiche                                   *
 *   $Id: strand.h,v 1.2 2008-01-24 10:20:10 wash Exp $           *
 *********************************************************************/


#ifndef STRAND_H
#define STRAND_H

#include <stdio.h>
#include <stdlib.h>

#define CODE_SCI 0
#define CODE_MEANMFE 1
#define CODE_CONSMFE 2
#define CODE_Z 3
#define CODE_SCI_Z 4
#define CODE_SCI_Z_MEANMFE 5
#define CODE_SCI_Z_CONSMFE 6
#define CODE_SCI_Z_MEANMFE_CONSMFE 7
#define CODE_MEANMFE_CONSMFE 8
#define CODE_Z_MEANMFE_CONSMFE 9
#define CODE_Z_CONSMFE 10
#define CODE_SCI_CONSMFE 11
#define CODE_SCI_MEANMFE 12
#define CODE_SCI_MEANMFE_CONSMFE 13
#define CODE_Z_MEANMFE 14

/* Get values for strand predictors */
void get_strand_predictors(
    char* window[], 
    char* window_reverse[], 
    char* names[],
    int length,
    int length_rev,
    int n_seq, 
    double *deltaSCI, 
    double *deltaMeanMFE, 
    double *deltaConsMFE,
    double *deltaZ,
    double *consFreqGU,
    char *modelDir);

/* Predicts reading direction of structured RNA candidate */
int strand(
    double deltaSCI, 
    double deltaMeanMFE, 
    double deltaConsMFE, 
    double deltaZ, 
    int n_seq, 
    double id, 
    double consFreqGU, 
    double* score,
    double* prob, 
    double *decValue, 
    char *modelDir,
    int code);

/* GU base pair frequency */
double frequency_gu_pairs(char *string, char *struc);

/* Consensus of alignment */
static char *strand_consensus(const char *AS[]);

/****************************************************************************************************************
  Read and scale descriptor values.
*****************************************************************************************************************/
static int read_scale_vector(FILE *scaleFile);
/* checks if the computet descriptors have their values in those of the decision.model 
	returns 1 if one or more of the descriptors is out of Range
			  0 otherwise 
	error message will be print into 'char * error' */
int checkStrandDescriptors(double* descriptors, int n_descr, char* error, int code);
void scale_strand_node_linear(struct svm_node* node,double scale[][3]);
void scale_strand_decision_node(struct svm_node* node);

/****************************************************************************************************************
  Init, read and free SVM model. 
*****************************************************************************************************************/
void strand_svm_init(char* basefilename);
void get_strand_model(char* basefilename);
static void strand_svm_free();

/****************************************************************************************************************
  Prediction 
*****************************************************************************************************************/
void predict_strand(double* prob, double* decValue, double* descriptors, int n_seq, double id, double consFreqGU, int n_descr);

/****************************************************************************************************************
  Functions to get prediction values for each strand
*****************************************************************************************************************/
double get_strand_forward_SCI();
double get_strand_forward_meanMFE();
double get_strand_forward_consMFE();
double get_strand_forward_Z();
double get_strand_reverse_SCI();
double get_strand_reverse_meanMFE();
double get_strand_reverse_consMFE();
double get_strand_reverse_Z();
char*  get_strand_output();
char*  get_strand_warning();

/****************************************************************************************************************
  Functions which free dynamic global variables.
*****************************************************************************************************************/
void free_strand();

/****************************************************************************************************************
  Miscellaneous
*****************************************************************************************************************/
void strand_env_variable_init(char* envVariable);
void code2string(int code);
void i2descriptor(int code, int i, char* descr);

#endif
