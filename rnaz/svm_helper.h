/*********************************************************************
 *                                                                   *
 *                          svm_helper.h                             *
 *                                                                   *
 *   Functions relating to SVM regression/classification and         *
 *    for interacting with the SVMLIB libraries                      *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                  *
 *                                                                   *
 *	   $Id: svm_helper.h,v 1.3 2004-09-19 13:31:41 wash Exp $    *
 *                                                                   *
 *********************************************************************/


/*
  declare svm_model for other functions (declaration is only in
  svmlib/svm.cpp and not in svmlib/svm.h)
*/

struct svm_model
{
  struct svm_parameter param;	// parameter 
  int nr_class; // number of classes, = 2 in regression/one class svm
  int l;			// total #SV
  struct svm_node **SV;		// SVs (SV[l])
  double **sv_coef;	// coefficients for SVs in decision functions (sv_coef[n-1][l])
  double *rho;		// constants in decision functions (rho[n*(n-1)/2])
  double *probA;          // pariwise probability information
  double *probB;

  // for classification only

  int *label;		// label of each class (label[n])
  int *nSV;		// number of SVs for each class (nSV[n])
				// nSV[0] + nSV[1] + ... + nSV[n-1] = l
  // XXX
  int free_sv;		// 1 if svm_model is created by svm_load_model
  // 0 if svm_model is created by svm_train
};


void get_regression_models(struct svm_model** avg_model,
						   struct svm_model** stdv_model,
						   char *basefilename);

struct svm_model* get_decision_model(char *basefilename);

struct svm_model* default_avg_model();

struct svm_model* default_stdv_model();

struct svm_model* default_decision_model();

void scale_regression_node(struct svm_node* node);

void scale_decision_node(struct svm_node* node);

void backscale_regression(double* avg, double* stdv);
