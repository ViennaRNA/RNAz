#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "svm.h"
#include "svm_helper.h"


/* All values in the SVMs are normalized to mean 0 and standard
   deviation 1. Therefore the mean and standard deviation of the
   training set must be known to scale the values appropriately.*/

/*                              Mean  Std.Dev. */
double scale_regression[][2]={{0.5,0.1581213081},   /* GC-content */
							  {0.5,0.1581213081},  /* "A-content" */
							  {0.5,0.1581213081},  /* "C-content" */
							  {225.0,114.569772373}};  /* length */


/* double scale_decision[][2]= {{-7.6,2.21 }, */
/* 							 {0,1.23 }, */
/* 							 {0,99.55},    */
/* 							 {0,6}};       */

double scale_decision[][2]= {{-7.87,2.76 },/* z */
							 {0,1.23 },   /* SCI */
							 {52.35,99.55},   /* ID */
							 {2,6}};      /* k  */


/* Also y-values are scaled in the regression SVMs, so also here the
   mean and standard deviation are needed */

double scale_avg[]={-58.60276,45.24618};
double scale_stdv[]={4.098457,1.107606};



/* Scales a SVM-lib "node" according to the scale_arrays defined
   above. Index 1 of the node is scaled with values of index 0 of the
   array.*/

void scale_node(struct svm_node* node,double scale[][2]){

  int i=0;
  double mu,sigma;
  while (1){
	mu=scale[node[i].index-1][0];
	sigma=scale[node[i].index-1][1];
	node[i].value=(node[i].value-mu)/sigma;
	++i;
	if (node[i].index==-1) break;
  }
}

void scale_node_linear(struct svm_node* node,double scale[][2]){

  int i=0;
  double min,max;
  double from,to;
  from=-1;
  to=1;
  while (1){
	min=scale[node[i].index-1][0];
	max=scale[node[i].index-1][1];
	node[i].value=from+(to-from)*(node[i].value-min)/(max-min);
	++i;
	if (node[i].index==-1) break;
  }
}



/* Scales x-data for the z-score regression */

void scale_regression_node(struct svm_node* node){
  scale_node(node,scale_regression);
}


void scale_decision_node(struct svm_node* node){
  scale_node_linear(node,scale_decision);
}



/* Back-scales y-data to original dimensions */

void backscale_regression(double* avg, double* stdv){
  *avg=(*avg)*scale_avg[1]+scale_avg[0];
  *stdv=(*stdv)*scale_stdv[1]+scale_stdv[0];
}


/* Loads both models for average and standard deviation from files
   given by a common basename. If no name is given, default models
   hard-coded in this file are used */

void get_regression_models(struct svm_model** avg_model, struct svm_model** stdv_model, char *basefilename){

  char fn[256];

  if (basefilename!=NULL){
	strncpy(fn, basefilename, 238); strcat(fn, "/mfe_avg.model");
	*avg_model = svm_load_model(fn);
	//if (*avg_model!=NULL) {
	//  printf("Loaded: %s (%p)\n",fn,*avg_model);
	//}
	strncpy(fn, basefilename, 237); strcat(fn, "/mfe_stdv.model");
	*stdv_model = svm_load_model(fn);
	//if (*stdv_model!=NULL) printf("Loaded: %s\n",fn);
  } else {
	*avg_model=default_avg_model();
	*stdv_model=default_stdv_model();
  }
 
}


struct svm_model* get_decision_model(char *basefilename){

  char fn[256];

  struct svm_model* model=NULL;
  
  if (basefilename!=NULL){
	strcpy(fn, basefilename); strcat(fn, "/decision.model");

	model=svm_load_model(fn);

  } else {
	model=default_decision_model();
  }

  return model;
  
}



/* Defines and returns default average model */
/* Not implemented yet */

struct svm_model* default_avg_model(){

  struct svm_model* model;

/*   struct svm_parameter par; */
  
/*   par.svm_type=4; */
/*   par.kernel_type=2; */
/*   par.gamma=0.5; */
/*   par.degree=0; */
/*   par.coef0=0; */

/*   model.param=par; */

/*   model.nr_class=2; */
/*   model.l=2; */

/*   model.sv_coef[1][0]=-5; */
/*   model.sv_coef[1][1]=-5; */

/*   model.rho[0]=26.1496; */

/*   model.probA=NULL; */
/*   model.probB=NULL; */
 
/*   model.label=NULL; */
/*   model.nSV=NULL; */
 
/*   model.free_sv=1; */
 
/*   model.SV[0][0].index=1; model.SV[0][0].value=-1; */
/*   model.SV[0][1].index=2; model.SV[0][1].value=-1; */
/*   model.SV[0][2].index=3; model.SV[0][2].value=-0.7; */
/*   model.SV[0][3].index=4; model.SV[0][3].value=0.1; */

/*   model.SV[1][0].index=1; model.SV[1][0].value=-1; */
/*   model.SV[1][1].index=2; model.SV[1][1].value=-1; */
/*   model.SV[1][2].index=3; model.SV[1][2].value=-0.7; */
/*   model.SV[1][3].index=4; model.SV[1][3].value=-0.7; */


  model=NULL;
  
  return model;

  
}


/* Defines and returns default standard deviation model */
/* Not implemented yet */

struct svm_model* default_stdv_model(){

  struct svm_model* model;
  model=NULL;
  return model;
 
}

struct svm_model* default_decision_model(){

  struct svm_model* model;
  model=NULL;
  return model;
 
}
