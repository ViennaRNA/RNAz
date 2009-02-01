/*********************************************************************
 *                                                                   *
 *                              zscore.c                             *
 *                                                                   *
 *	Compute a z-score to assess significance of a predicted MFE      *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                      *
 *                                                                   *
 *	   $Id: zscore.c,v 1.6 2006/10/12 13:19:01 wash Exp $            *
 *                                                                   *
 *********************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "fold.h"
#include "svm.h"
#include "svm_helper.h"
#include "zscore.h"


struct svm_model *avg_model, *stdv_model;


/* Initializes pointers to the two regression models. If a basename is
   given, models are loaded from files, otherwise standard models are
   used (not implemented yet) */

void regression_svm_init(char *basefilename){

  avg_model=NULL;
  stdv_model=NULL;

  get_regression_models(&avg_model,&stdv_model,basefilename);

  if (avg_model==NULL)
	nrerror("ERROR: Could not load mu-regression model. \n");
  
  if (stdv_model==NULL)
	nrerror("ERROR: Could not load sigma-regression model.\n");
}

/* Destroys SVM-models */

void regression_svm_free(){

  svm_destroy_model(avg_model);
  svm_destroy_model(stdv_model);

  avg_model=NULL;
  stdv_model=NULL;
}


/* Predict mean and standard deviation of random sequences for a base
   composition of a given sequence seq*/

void predict_values(const char *seq, double *avg, double *stdv) {

  double GC,A,C;
  struct svm_node node[5];
  int i, length,n_A,n_C,n_T,n_G;

  n_A=n_T=n_G=n_C=0;
  
  for (i=0; i<strlen(seq); i++) {
    switch (seq[i]) {
    case 'A': n_A++; break;
    case 'C': n_C++; break;
    case 'G': n_G++; break;
	case 'T': n_T++; break;
	case 'U': n_T++; break;  
    }
  }

  length=strlen(seq);

  if (length==0 || n_A+n_T==0 || n_G+n_C==0){
	GC=A=C=0;
  } else {
	GC=((double)(n_G+n_C))/length;
	A=((double)n_A)/(n_A+n_T);
	C=((double)n_C)/(n_G+n_C);
  }

  node[0].index = 1; node[0].value = GC;
  node[1].index = 2; node[1].value = A;
  node[2].index = 3; node[2].value = C;
  node[3].index = 4; node[3].value = length;
  node[4].index =-1;

  scale_regression_node((struct svm_node*)&node);

  *avg=svm_predict(avg_model,node);
  *stdv=svm_predict(stdv_model,node);

  backscale_regression(avg,stdv);

}


/* Predict the z-score of a sequence. If a mfe>0 is given, the mfe of
   the sequence is calculated by fold() otherwise the precalculated
   value is used*/

double mfe_zscore(const char *seq, double mfe) {
  double E, stdv, avg;
  char *struc;

  if (mfe>0){
  
	struc = space(strlen(seq)+1);
	E = fold(seq, struc);
	free(struc);
  } else {

	E=mfe;
	
  }
  
  predict_values(seq, &avg, &stdv);

  return ((E-avg)/stdv);

}
  
