/*********************************************************************
 *                                                                   *
 *                          svm_helper.h                             *
 *                                                                   *
 *   Functions relating to SVM regression/classification and         *
 *    for interacting with the SVMLIB libraries                      *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                  *
 *                                                                   *
 *	   $Id: svm_helper.h,v 1.5 2006/01/29 18:22:37 wash Exp $    *
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

struct svm_model *svm_load_model_string(char *fp);

int print_model(const char *model_file_name, struct svm_model *model);



