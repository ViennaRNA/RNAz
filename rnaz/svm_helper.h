/*********************************************************************
 *                                                                   *
 *                          svm_helper.h                             *
 *                                                                   *
 *   Functions relating to SVM regression/classification and         *
 *    for interacting with the SVMLIB libraries                      *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                  *
 *                                                                   *
 *	   $Id: svm_helper.h,v 1.4 2004-09-22 10:05:21 wash Exp $    *
 *                                                                   *
 *********************************************************************/


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
