/*********************************************************************
 *                                                                   *
 *                              zscore.h                             *
 *                                                                   *
 *	Compute a z-score to assess significance of a predicted MFE      *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                      *
 *                                                                   *
 *	   $Id: zscore.h,v 1.2 2004-09-19 12:35:54 wash Exp $          *
 *                                                                   *
 *********************************************************************/

void regression_svm_init(char *basefilename);

void regression_svm_free();

double mfe_zscore(const char *seq,double mfe);
