#ifndef STRAND_H
#define STRAND_H

#include <stdio.h>
#include <stdlib.h>

/* Get values for strand predictors */
void get_strand_predictors(
    char* window[], 
    char* window_reverse[], 
    int length, 
    int n_seq, 
    double *deltaSCI, 
    double *deltaMeanMFE, 
    double *deltaConsMFE,
    double *deltaZ,
    char *modelDir);

/* Predicts reading direction of structured RNA candidate */
int predict_strand(
    double deltaSCI, 
    double deltaMeanMFE, 
    double deltaConsMFE, 
    double deltaZ, 
    int n_seq, 
    double id, 
    int* strand, 
    double* prob, 
    double *decValue, 
    char *modelDir);

/* checks if the computet descriptors have their values in those of the decision.model 
	returns 1 if one or more of the descriptors is out of Range
			  0 otherwise 
	error message will be print into 'char * error' */
int checkStrandDescriptors(double deltaSCI, double deltaMeanMFE, double deltaConsMFE, double deltaZ, char* error);
			
/* svm classification depending on various variables */
void predict_strand_SCI_MeanMFE_ConsMFE(double* prob, double* decValue, double deltaSCI, double deltaMeanMFE, double deltaConsMFE, int, double);
void predict_strand_SCI_MeanMFE_ConsMFE_Z(double* prob, double* decValue, double deltaSCI, double deltaMeanMFE, double deltaConsMFE, double delta, int, double);
void predict_strand_SCI_MeanMFE(double* prob, double* decValue, double deltaSCI, double deltaMeanMFE, int, double);
void predict_strand_SCI_ConsMFE(double* prob, double* decValue, double deltaSCI, double deltaConsMFE, int, double);
void predict_strand_SCI_Z(double* prob, double* decValue, double deltaSCI, double deltaZ, int, double);
void predict_strand_ConsMFE_Z(double* prob, double* decValue, double deltaConsMFE, double deltaZ, int n_seq, double id);
void predict_strand_SCI_ConsMFE_Z(double *prob, double *decValue, double deltaSCI, double deltaConsMFE, double deltaZ, int n_seq, double id);
void predict_strand_MEANMFE_CONSMFE(double *prob, double *decValue, double deltaMeanMFE, double deltaConsMFE, int n_seq, double id);
void predict_strand_MEANMFE_CONSMFE_Z(double *prob, double *decValue, double deltaMeanMFE, double deltaConsMFE, double deltaZ, int n_seq, double id);

void strand_env_variable_init(char* envVariable);
void strand_svm_init(char* basefilename);
void get_strand_model(char* basefilename);
void strand_svm_free();

void scale_strand_node_linear(struct svm_node* node,double scale[][3]);
void scale_strand_decision_node(struct svm_node* node);

#endif
