/*********************************************************************
 *                                                                   *
 *                              zscore.c                             *
 *                                                                   *
 *	Compute a z-score to assess significance of a predicted MFE  *
 *                                                                   *
 *	          c Stefan Washietl, Ivo L Hofacker                  *
 *                                                                   *
 *	   $Id: zscore.c,v 1.6 2006/10/12 13:19:01 wash Exp $        *
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
#include "fold_vars.h"

#define IN_RANGE(LOWER,VALUE,UPPER) ((VALUE <= UPPER) && (VALUE >= LOWER))

struct svm_model *avg_model, *stdv_model;
struct svm_model *GC20_30_avg, *GC30_36_avg, *GC36_40_avg, 
  *GC40_46_avg, *GC46_50_avg, *GC50_56_avg, *GC56_60_avg,
  *GC60_66_avg, *GC66_70_avg, *GC70_80_avg;
struct svm_model *GC20_30_stdv, *GC30_36_stdv, *GC36_40_stdv, 
  *GC40_46_stdv, *GC46_50_stdv, *GC50_56_stdv, *GC56_60_stdv,
  *GC60_66_stdv, *GC66_70_stdv, *GC70_80_stdv;


/* Calculates the base frequencies for both mononucleotides and 
   dinucleotides */

void base_frequencies(const char *array, unsigned int n, double *mono, double *di)
{
  unsigned int i;
  char dinuc[3];
  dinuc[2] = '\0';
  
  for (i = 0; i < n ; i++)
  {
    switch( array[i] )
    {
    case 'A' : mono[0]++;break;
    case 'C' : mono[1]++;break;
    case 'G' : mono[2]++;break;
    case 'T' : mono[3]++;break;
    case 'U' : mono[3]++;break; 
    default  : mono[4]++;break;
    }
    dinuc[0] = array[i];
    dinuc[1] = array[i+1];
    
    
    if (strcmp(dinuc,"AA")==0) di[0]++;
    if (strcmp(dinuc,"AC")==0) di[1]++;
    if (strcmp(dinuc,"AG")==0) di[2]++;
    if (strcmp(dinuc,"AT")==0) di[3]++;
    if (strcmp(dinuc,"AU")==0) di[3]++;
    if (strcmp(dinuc,"CA")==0) di[4]++;
    if (strcmp(dinuc,"CC")==0) di[5]++;
    if (strcmp(dinuc,"CG")==0) di[6]++;
    if (strcmp(dinuc,"CT")==0) di[7]++;
    if (strcmp(dinuc,"CU")==0) di[7]++;
    if (strcmp(dinuc,"GA")==0) di[8]++;
    if (strcmp(dinuc,"GC")==0) di[9]++;
    if (strcmp(dinuc,"GG")==0) di[10]++;
    if (strcmp(dinuc,"GT")==0) di[11]++;
    if (strcmp(dinuc,"GU")==0) di[11]++;
    if (strcmp(dinuc,"TA")==0) di[12]++;
    if (strcmp(dinuc,"TC")==0) di[13]++;
    if (strcmp(dinuc,"TG")==0) di[14]++;
    if (strcmp(dinuc,"TT")==0) di[15]++;
    if (strcmp(dinuc,"UA")==0) di[12]++;
    if (strcmp(dinuc,"UC")==0) di[13]++;
    if (strcmp(dinuc,"UG")==0) di[14]++;
    if (strcmp(dinuc,"UU")==0) di[15]++;
  }

  /* the last one */
  switch( array[n] )
  {
    case 'A' : mono[0]++;break;
    case 'C' : mono[1]++;break;
    case 'G' : mono[2]++;break;
    case 'T' : mono[3]++;break;
    case 'U' : mono[3]++;break; 
    default  : mono[4]++;break;
  }
  
  for (i = 0; i < 5 ; i++)
  {
    char tmp[10];
    /*sprintf(tmp, "%.3f", (double) mono[i]/n);
      mono[i] = (double) atof(tmp);*/
    mono[i] = (double) mono[i]/n;
  }

  for (i = 0; i < 16 ; i++)
  {
    char tmp[10];
    sprintf(tmp, "%.3f", (double) di[i]/(n-1));
    di[i] = (double) atof(tmp);
  }

}

/* Uses a Fisher-Yates shuffle to generate a mononucleotide
   shuffled sequence. */

void fisher_yates_shuffle(char *array, size_t n, char *work_array)
{
  size_t k;
  for (k = 0; k < n ; k++) {
    work_array[k] = array[k];
  }
  
  if (n > 1) {
    size_t i;
    for (i = 0; i < n - 1; i++) {
      size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}

/* Uses a Altschul-Erickson shuffle to generate a
   dinucleotide shuffled sequence. */
void altschul_erickson_shuffle(char *array, size_t n, char *work_array)
{
  int **E_s;
  int **Z;
  int *n_E;
  int *i_E;
  size_t i,j;
  int pre;
  int eulerian_flag = 0;
  int s_1;
  int s_f;

  /* transform string array to numbers and initialize some arrays*/
  for (i = 0; i < n ; i++)
  {
    switch( array[i] )
    {
    case 'A' : work_array[i] = 0;break;
    case 'C' : work_array[i] = 1;break;
    case 'G' : work_array[i] = 2;break;
    case 'T' : work_array[i] = 3;break;
    case 'U' : work_array[i] = 3;break; 
    default  : work_array[i] = 4;break;
    }
  }
  s_1 = work_array[0];
  s_f = work_array[n-1];

  Z  = (int **) space(sizeof(int *) * 5);
  for (i = 0; i < 5; i++)
  {
    Z[i]  = (int *) space(sizeof(int) * 5);
  }
  

  /* (1) Construct the doublet graph G and the edge ordering E corresponding to S.*/
  E_s  = (int **) space(sizeof(int *) * 5);
  n_E = (int *) space(sizeof(int)    * 5);
  i_E = (int *) space(sizeof(int)    * 5);
  for (i = 0; i < 5; i++)
  {
    E_s[i]  = (int *) space(sizeof(int) * (n-1));
    n_E[i] = 0;
    i_E[i] = 0;
  }
  

  pre = (int) work_array[0];
  for (i = 1; i < n; i++)
  {
    E_s[pre][n_E[pre]] = (int) work_array[i];
    n_E[pre]++;
    pre = (int) work_array[i];
  }
  /* printf("A:%d, C:%d, G:%d, T:%d, N:%d\n",n_E[0],n_E[1],n_E[2],n_E[3],n_E[4]);*/

  while (eulerian_flag == 0)
  {
    /* (2) For each vertex s in G except s_f, randomly select one edge from the s
       edge list of E(S) to be the last edge of the s list in a new edge ordering.*/
    for (i = 0; i < 5; i++)
    {
      int rand_index;
      int swap;
      /*printf("i:%d, #:%d\n",i,n_E[i]);*/
      if (n_E[i] == 0 || i == s_f) continue;
      rand_index = (int)(rand() / (((double)RAND_MAX + 1)/ n_E[i]));
      /*swap*/
      swap = E_s[i][n_E[i]-1];
      E_s[i][n_E[i]-1] = E_s[i][rand_index];
      E_s[i][rand_index] = swap;
    }
    

    /* (3) From this set of last edges, construct the last-edge graph Z and
       determine wheter or not all of its vertices are connected to s_f.*/

    /*set everything to zero in Z*/
    for (i = 0; i < 5; i++)
    {
      for (j = 0; j < 5; j++)
      {
	Z[i][j] = 0;
      }
    }
    

    int max = 0;
    for (i = 0; i < 5; i++)
    {
      if (n_E[i] > 0)
      {
	Z[i][E_s[i][n_E[i]-1]] = 1;
	max++;
      }
    }
    
    
    if (n_E[s_f] == 0)
    {
      Z[s_f][s_f] = 1;
      max++;
    }

    int seen = 0;
    for (i = 0; i < 5; i++)
    {
      int node = i;
      int dummy;
      for (dummy = 0; dummy < max; dummy++)
      {
	int connected = -1;
	for (j = 0; j < 5; j++)
	{
	  if (Z[node][j] == 1) connected = j;
	}
	
	/*node is not connected to any other node, we can stop*/
	if (connected == -1)
	{
	  break;
	}
	/*is it s_f?*/
	if (connected == s_f)
	{
	  seen++;
	  break;
	}
	/*now search for the connecred node if it is connected to s_f*/
	node = connected;
      }
    }
    if (seen == max) eulerian_flag = 1;
    
    
    /* (4) If any vertex is not connected in Z to s_f, the new edge ordering
       will not be eulerian, so return to (2). If all vertices are connected
       in Z to s_f, the new edge ordering will be Eulerian, so continue to (5)*/
  }

  /* (5) For each vertex s in G, randomly permute the remaining edges of the s
     edge list of E(S) to generate the s edge list of the new edge ordering E(S').*/
  for (i = 0; i < 5; i++)
  {
    if (n_E[i] > 0)
    {
	for (j = 0; j < n_E[i]-1; j++)
	{
	  int rand_index = (int)(rand() / (((double)RAND_MAX + 1)/ (n_E[i]-1)));
	  int tmp = E_s[i][j];
	  E_s[i][j] = E_s[i][rand_index];
	  E_s[i][rand_index] = tmp;
	}
    }
  }
  
  
  /*(6) Construct sequence S', a random DP permutation of S, from E(S') as follows.
    Start at the     s_1 edge list. At each s_i edge list, add s_i to S', delete the
    first edge s_is_j of the edge list, and move to the s_j edge list. Continue this
    process until all edge lists are exhausted.*/

  i = 0;
  j = s_1;
  int flag = 0;
  while (flag == 0)
  {
    int tmp;
    switch( j )
    {
    case 0 : work_array[i] = 'A';break;
    case 1 : work_array[i] = 'C';break;
    case 2 : work_array[i] = 'G';break;
    case 3 : work_array[i] = 'U';break;
    case 4 : work_array[i] = 'N';break;  
    }
    tmp = E_s[j][i_E[j]];
    i_E[j]++;
    j = tmp;
    i++;
    if (i_E[j] == n_E[j]) flag = 1;
  }
  /*now set the last*/
  switch( s_f )
  {
    case 0 : work_array[i] = 'A';break;
    case 1 : work_array[i] = 'C';break;
    case 2 : work_array[i] = 'G';break;
    case 3 : work_array[i] = 'U';break;
    case 4 : work_array[i] = 'N';break;  
  }
  

  /* free all the stuff */
  for (i = 0; i < 5; i++)
  {
    free(E_s[i]);
    free(Z[i]);
  }

  free(E_s);
  free(Z);
  free(n_E);
  free(i_E);
}

void zscore_explicitly_shuffled(const char *seq, double *avg, double *stdv, int type){
    unsigned int n = 1000;
    unsigned int counter;
    unsigned int length = strlen(seq);
    char *tmp = (char *) space((unsigned) length+1);
    char *shuff = (char *) space((unsigned) length+1);
    char *structure = (char *) space((unsigned) length+1);
    float *energies = (float*) space(sizeof(float) * n);
    float mean = 0;
    float sd = 0;
    /* generate seed from time */
    srand((unsigned)time(NULL));
    strcpy(tmp, seq);

    /*printf("INHERE\n");*/
    
    for (counter = 0; counter < n; counter++)
    {
      if (type == 1) fisher_yates_shuffle(tmp,length,shuff);
      if (type == 3) altschul_erickson_shuffle(tmp,length,shuff);
      energies[counter] = fold(shuff, structure);
      free_arrays();
      mean = mean + energies[counter];
    }

    mean = (float) mean/n;
    
    for (counter = 0; counter < n; counter++)
    {
      sd += (energies[counter]-mean)*(energies[counter]-mean);
    }

    sd = (float) sd/(n-1);
    sd = sqrt(sd);

    *avg = (double) mean;
    *stdv = (double) sd;

    free(energies);
    free(tmp);
    free(shuff);
    free(structure);
}



/* Initializes pointers to the two regression models. If a basename is
   given, models are loaded from files, otherwise standard models are
   used (not implemented yet) */

void regression_svm_init(){

  avg_model=NULL;
  stdv_model=NULL;

  get_regression_models(&avg_model,&stdv_model,-1);

  if (avg_model==NULL)
  	nrerror("ERROR: Could not load mu-regression model. \n");
  
  if (stdv_model==NULL)
	nrerror("ERROR: Could not load sigma-regression model.\n");
}

/* Destroys SVM-models */

void regression_svm_free(){

  if (avg_model != NULL) svm_destroy_model(avg_model);
  if (stdv_model != NULL) svm_destroy_model(stdv_model);
  if (GC20_30_avg != NULL) svm_destroy_model(GC20_30_avg);
  if (GC30_36_avg != NULL) svm_destroy_model(GC30_36_avg);
  if (GC36_40_avg != NULL) svm_destroy_model(GC36_40_avg);
  if (GC40_46_avg != NULL) svm_destroy_model(GC40_46_avg);
  if (GC46_50_avg != NULL) svm_destroy_model(GC46_50_avg);
  if (GC50_56_avg != NULL) svm_destroy_model(GC50_56_avg);
  if (GC56_60_avg != NULL) svm_destroy_model(GC56_60_avg);
  if (GC60_66_avg != NULL) svm_destroy_model(GC60_66_avg);
  if (GC66_70_avg != NULL) svm_destroy_model(GC66_70_avg);
  if (GC70_80_avg != NULL) svm_destroy_model(GC70_80_avg);
  if (GC20_30_stdv != NULL) svm_destroy_model(GC20_30_stdv);
  if (GC30_36_stdv != NULL) svm_destroy_model(GC30_36_stdv);
  if (GC36_40_stdv != NULL) svm_destroy_model(GC36_40_stdv);
  if (GC40_46_stdv != NULL) svm_destroy_model(GC40_46_stdv);
  if (GC46_50_stdv != NULL) svm_destroy_model(GC46_50_stdv);
  if (GC50_56_stdv != NULL) svm_destroy_model(GC50_56_stdv);
  if (GC56_60_stdv != NULL) svm_destroy_model(GC56_60_stdv);
  if (GC60_66_stdv != NULL) svm_destroy_model(GC60_66_stdv);
  if (GC66_70_stdv != NULL) svm_destroy_model(GC66_70_stdv);
  if (GC70_80_stdv != NULL) svm_destroy_model(GC70_80_stdv);

  avg_model=NULL;
  stdv_model=NULL;
}


/* Predict mean and standard deviation of random sequences for a base
   composition of a given sequence seq*/
/* type = 0: use MONO-nucleotide shuffled SVM */
/* type = 1: explictily shuffle MONO-nucleotide */
/* type = 2: use DI-nucleotide shuffled SVM */
/* type = 3: explictily shuffle DI-nucleotide */

void predict_values(const char *seq, double *avg, double *stdv, int *type, 
		    int avoid_shuffle, char* warning_string) {

  unsigned int counter;
  double *mono_array = (double*) space(sizeof(double) * 5);
  double *di_array = (double*) space(sizeof(double) * 16);
  unsigned int length = strlen(seq);
  double GplusC,AT_ratio,CG_ratio;
  char tmp[10];
  int verbose;
  verbose = 1; /* set to 1 to get out of range warnings. */

  /* count base frequencies */
  /* was allocated with space, which makes a calloc */
  for (counter = 0; counter < 5; counter++) { mono_array[counter] = 0; }
  for (counter = 0; counter < 16; counter++) { di_array[counter] = 0; }
  base_frequencies(seq, length, mono_array, di_array);

  /* RNAz 1.0 uses always the full float. In RNAz 2.0 we used for training numbers 
     with at most 3 decimal places. We adujst it here for RNAz 2.0. */
  if (*type > 0) {
    for (counter = 0; counter < 5; counter++) { 
      char tmp[10];
      sprintf(tmp, "%.3f", (double) mono_array[counter]);
      mono_array[counter] = (double) atof(tmp);
    }
  }

  /* calculate descriptors required both for mono and dinucleotide regression */
  /* G+C content */

  if (*type == 0) {
    GplusC =(double) mono_array[1] + mono_array[2];
  } else {
    /* Once again we round to three decimal places for RNAz 2.0 */
    sprintf(tmp, "%.3f", (double) mono_array[1] + mono_array[2]);
    GplusC = (double) atof(tmp);
  }
  if (mono_array[0] + mono_array[3]==0 || mono_array[1] + mono_array[2]==0) {
    AT_ratio = CG_ratio = 0.0;
  } else {
    if (*type == 0) {
      CG_ratio = (double) mono_array[1] / (mono_array[1] + mono_array[2]);
      AT_ratio = (double) mono_array[0] / (mono_array[0] + mono_array[3]);
    } else {
      /* Once again we round to three decimal places for RNAz 2.0 */
      /* CG ratio = C / (C + G) */
      sprintf(tmp, "%.3f", (double) mono_array[1] / (mono_array[1] + mono_array[2]));
      CG_ratio = (double) atof(tmp);
      /* AT ratio = A / (A + T) */
      sprintf(tmp, "%.3f", (double) mono_array[0] / (mono_array[0] + mono_array[3]));
      AT_ratio = (double) atof(tmp);
    }
  }

  /* let's check the bounds so we can decide what to do */
  if (*type==0) {
    if (dangles != 2) {
      if (verbose == 1) {
	strcpy(warning_string," WARNING: Regression was not trained with dangles option other than -d2.\n");
	warning_string+=strlen(warning_string);
      }
      *type = 1;
    }
    if (length < 50) {
      if (verbose == 1) {
	strcpy(warning_string," WARNING: Sequence too short.\n");
	warning_string+=strlen(warning_string);
      }
      *type = 1;
    }
    if (length > 400) {
      if (verbose == 1) {
	strcpy(warning_string," WARNING: Sequence too long.\n");
	warning_string+=strlen(warning_string);
      }
      *type = 1;
    }
    if ((!IN_RANGE(0.25,GplusC,0.75)) ||
	(!IN_RANGE(0.25,AT_ratio,0.75)) ||
	(!IN_RANGE(0.25,CG_ratio,0.75))) {
      if (verbose == 1) {
	strcpy(warning_string," WARNING: Base composition out of range.\n");
	warning_string+=strlen(warning_string);
      }
      *type = 1;
    }
  }
  
  if (*type==2) {
    if (dangles != 2) {
      if (verbose == 1) {
	strcpy(warning_string," WARNING: Regression was not trained with dangles option other than -d2.\n");
	warning_string+=strlen(warning_string);
      }
      *type = 3;
    }
    if (length < 50) {
      if (verbose == 1) {
	strcpy(warning_string," WARNING: Sequence too short.\n");
	warning_string+=strlen(warning_string);
      }
      *type = 3;
    }
    if (length > 200) {
      if (verbose == 1) {
	strcpy(warning_string," WARNING: Sequence too long.\n");
	warning_string+=strlen(warning_string);
      }
      *type = 3;
    }

    if ((!IN_RANGE(0.20,GplusC,0.80)) ||
	(!IN_RANGE(0.20,AT_ratio,0.80)) ||
	(!IN_RANGE(0.20,CG_ratio,0.80))) {
      if (verbose == 1) {
	strcpy(warning_string," WARNING: Base composition out of range.\n");
	warning_string+=strlen(warning_string);
	strcpy(warning_string," WARNING: Base composition out of range.\n");
	warning_string+=strlen(warning_string);
      }
       *type = 3;
    }

    /* let's see if we are in bounds for dinculeotide params */
    /* AA */
    if (di_array[0] > 0) {
      if(abs((di_array[0]-mono_array[0]*mono_array[0])/(-mono_array[0]*mono_array[0])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* AC */
    if (di_array[1] > 0) {
      if(abs((di_array[1]-mono_array[0]*mono_array[1])/(-mono_array[0]*mono_array[1])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* AG */
    if (di_array[2] > 0) {
      if(abs((di_array[2]-mono_array[0]*mono_array[2])/(-mono_array[0]*mono_array[2])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* AT */
    if (di_array[3] > 0) {
      if(abs((di_array[3]-mono_array[0]*mono_array[3])/(-mono_array[0]*mono_array[3])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* CA */
    if (di_array[4] > 0) {
      if(abs((di_array[4]-mono_array[1]*mono_array[0])/(-mono_array[1]*mono_array[0])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* CC */
    if (di_array[5] > 0) {
      if(abs((di_array[5]-mono_array[1]*mono_array[1])/(-mono_array[1]*mono_array[1])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* CG */
    if (di_array[6] > 0) {
      if(abs((di_array[6]-mono_array[1]*mono_array[2])/(-mono_array[1]*mono_array[2])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* CT */
    if (di_array[7] > 0) {
      if(abs((di_array[7]-mono_array[1]*mono_array[3])/(-mono_array[1]*mono_array[3])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* GA */
    if (di_array[8] > 0) {
      if(abs((di_array[8]-mono_array[2]*mono_array[0])/(-mono_array[2]*mono_array[0])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* GC */
    if (di_array[9] > 0) {
      if(abs((di_array[9]-mono_array[2]*mono_array[1])/(-mono_array[2]*mono_array[1])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* GG */
    if (di_array[10] > 0) {
      if(abs((di_array[10]-mono_array[2]*mono_array[2])/(-mono_array[2]*mono_array[2])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* GT */
    if (di_array[11] > 0) {
      if(abs((di_array[11]-mono_array[2]*mono_array[3])/(-mono_array[2]*mono_array[3])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    } 
    /* TA */
    if (di_array[12] > 0) {
      if(abs((di_array[12]-mono_array[3]*mono_array[0])/(-mono_array[3]*mono_array[0])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* TC */
    if (di_array[13] > 0) {
      if(abs((di_array[13]-mono_array[3]*mono_array[1])/(-mono_array[3]*mono_array[1])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* TG */
    if (di_array[14] > 0) {
      if(abs((di_array[14]-mono_array[3]*mono_array[2])/(-mono_array[3]*mono_array[2])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    }
    /* TT */
    if (di_array[15] > 0) {
      if(abs((di_array[15]-mono_array[3]*mono_array[3])/(-mono_array[3]*mono_array[3])) > 1.5) {
	if (verbose == 1) {
	  strcpy(warning_string," WARNING: Dinucleotide composition out of range.\n");
	  warning_string+=strlen(warning_string);
	}
	*type = 3;
      }
    } 
  }

  /*printf("Type: %d\n", *type);*/
  /* now we have to see if user allowed explicite shuffling */
  if (*type == 1 && (length >= 50 && length <= 400 && avoid_shuffle)) *type = 0; /* mono */
  if (*type == 3 && (length >= 50 && length <= 400 && avoid_shuffle)) *type = 2; /* di */

  /*********************************************/
  /* now we DEFINITELY know which type to take */

  /*printf("Type: %d,%f,%f,%f,%f,%f,%f,%f\n", *type,mono_array[0],mono_array[1],mono_array[2],mono_array[3],GplusC,AT_ratio,CG_ratio);*/

  /* Mononucleotide Regression */
  if (*type == 0) {

    struct svm_node node_mono[5];

    /*printf("GplusC:%f,AT_ratio:%f,CG_ratio:%f,%d\n",GplusC,AT_ratio,CG_ratio,length);*/

    node_mono[0].index = 1; node_mono[0].value = GplusC;
    node_mono[1].index = 2; node_mono[1].value = AT_ratio;
    node_mono[2].index = 3; node_mono[2].value = CG_ratio;
    node_mono[3].index = 4; node_mono[3].value = length;
    node_mono[4].index =-1;
    
    scale_regression_node((struct svm_node*)&node_mono);
    
    *avg=svm_predict(avg_model,node_mono);
    *stdv=svm_predict(stdv_model,node_mono);

    backscale_regression(avg,stdv);
  }

  /* Mononucleotide explictly shuffled */
  if (*type == 1) {

    zscore_explicitly_shuffled(seq, avg, stdv, *type);
  }

  /* Dinucleotide Regression */
  if (*type == 2 ) {  

    double norm_length;
    struct svm_node node_di[21];
    
    /* normalized, scaled sequence length */
    sprintf(tmp, "%.5f", (double) (length-50)/150);
    norm_length = (double) atof(tmp);
    
    /* now set node properties */
    node_di[0].index = 1; node_di[0].value = GplusC;
    node_di[1].index = 2; node_di[1].value = CG_ratio;
    node_di[2].index = 3; node_di[2].value = AT_ratio;
    node_di[3].index = 4; node_di[3].value = di_array[0];
    node_di[4].index = 5; node_di[4].value = di_array[1];
    node_di[5].index = 6; node_di[5].value = di_array[2];
    node_di[6].index = 7; node_di[6].value = di_array[3];
    node_di[7].index = 8; node_di[7].value = di_array[4];
    node_di[8].index = 9; node_di[8].value = di_array[5];
    node_di[9].index = 10; node_di[9].value = di_array[6];
    node_di[10].index = 11; node_di[10].value = di_array[7];
    node_di[11].index = 12; node_di[11].value = di_array[8];
    node_di[12].index = 13; node_di[12].value = di_array[9];
    node_di[13].index = 14; node_di[13].value = di_array[10];
    node_di[14].index = 15; node_di[14].value = di_array[11];
    node_di[15].index = 16; node_di[15].value = di_array[12];
    node_di[16].index = 17; node_di[16].value = di_array[13];
    node_di[17].index = 18; node_di[17].value = di_array[14];
    node_di[18].index = 19; node_di[18].value = di_array[15];
    node_di[19].index = 20; node_di[19].value = norm_length;
    node_di[20].index =-1;

    /* Now we have to fetch the right model and calculate avg and stdv*/
  
    if (GplusC >= 0.200 && GplusC < 0.300) {
      if (GC20_30_avg == NULL || GC20_30_stdv == NULL)
      {
	get_regression_models(&GC20_30_avg, &GC20_30_stdv, 0);
	/*printf("\n==Model for GC 20-30 was loaded.\n");*/
      }
      *avg=svm_predict(GC20_30_avg,node_di);
      *stdv=svm_predict(GC20_30_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.300 && GplusC < 0.360) {
      if (GC30_36_avg == NULL || GC30_36_stdv == NULL)
      {
	get_regression_models(&GC30_36_avg, &GC30_36_stdv, 1);
	/*printf("\n==Model for GC 30-36 was loaded.\n");*/
      }
      *avg=svm_predict(GC30_36_avg,node_di);
      *stdv=svm_predict(GC30_36_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.360 && GplusC < 0.400) {
      if (GC36_40_avg == NULL || GC36_40_stdv == NULL)
      {
	get_regression_models(&GC36_40_avg, &GC36_40_stdv, 2);
	/*printf("\n==Model for GC 36-40 was loaded.\n");*/
      }
      *avg=svm_predict(GC36_40_avg,node_di);
      *stdv=svm_predict(GC36_40_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.400 && GplusC < 0.460) {
      if (GC40_46_avg == NULL || GC40_46_stdv == NULL)
      {
	get_regression_models(&GC40_46_avg, &GC40_46_stdv, 3);
	/*printf("\n==Model for GC 40-46 was loaded.\n");*/
      }
      *avg=svm_predict(GC40_46_avg,node_di);
      *stdv=svm_predict(GC40_46_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.460 && GplusC < 0.500) {
      if (GC46_50_avg == NULL || GC46_50_stdv == NULL)
      {
	get_regression_models(&GC46_50_avg, &GC46_50_stdv, 4);
	/*printf("\n==Model for GC 46-50 was loaded.\n");*/
      }
      *avg=svm_predict(GC46_50_avg,node_di);
      *stdv=svm_predict(GC46_50_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.500 && GplusC < 0.560) {
      if (GC50_56_avg == NULL || GC50_56_stdv == NULL)
      {
	get_regression_models(&GC50_56_avg, &GC50_56_stdv, 5);
	/*printf("\n==Model for GC 50-56 was loaded.\n");*/
      }
      *avg=svm_predict(GC50_56_avg,node_di);
      *stdv=svm_predict(GC50_56_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.560 && GplusC < 0.600) {
      if (GC56_60_avg == NULL || GC56_60_stdv == NULL)
      {
	get_regression_models(&GC56_60_avg, &GC56_60_stdv, 6);
	/*printf("\n==Model for GC 56-60 was loaded.\n");*/
      }
      *avg=svm_predict(GC56_60_avg,node_di);
      *stdv=svm_predict(GC56_60_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.600 && GplusC < 0.660) {
      if (GC60_66_avg == NULL || GC60_66_stdv == NULL)
      {
	get_regression_models(&GC60_66_avg, &GC60_66_stdv, 7);
	/*printf("\n==Model for GC 60-66 was loaded.\n");*/
      }
      *avg=svm_predict(GC60_66_avg,node_di);
      *stdv=svm_predict(GC60_66_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.660 && GplusC < 0.700) {
      if (GC66_70_avg == NULL || GC66_70_stdv == NULL)
      {
	get_regression_models(&GC66_70_avg, &GC66_70_stdv, 8);
	/*printf("\n==Model for GC 66-70 was loaded.\n");*/
      }
      *avg=svm_predict(GC66_70_avg,node_di);
      *stdv=svm_predict(GC66_70_stdv,node_di);
      *avg = *avg/10.0 * length;
    }

    if (GplusC >= 0.700 && GplusC <= 0.800) {
      if (GC70_80_avg == NULL || GC70_80_stdv == NULL)
      {
	get_regression_models(&GC70_80_avg, &GC70_80_stdv, 9);
	/*printf("\n==Model for GC 70-80 was loaded.\n");*/
      }
      *avg=svm_predict(GC70_80_avg,node_di);
      *stdv=svm_predict(GC70_80_stdv,node_di);
      *avg = *avg/10.0 * length;
    }
  }

  /* Dinucleotide explictly shuffled */
  if (*type == 3) {
    zscore_explicitly_shuffled(seq, avg, stdv, *type);
  }
  
  free(mono_array);
  free(di_array);
  
}


/* Predict the z-score of a sequence. If a mfe>0 is given, the mfe of
   the sequence is calculated by fold() otherwise the precalculated
   value is used*/
/* type = 0: use MONO-nucleotide shuffled SVM */
/* type = 1: explictily shuffle MONO-nucleotide */
/* type = 2: use DI-nucleotide shuffled SVM */
/* type = 3: explictily shuffle DI-nucleotide */

double mfe_zscore(const char *seq, double mfe, int *type, int avoid_shuffle,
		  char* warning_string) {
  double E, stdv, avg;
  char *struc;

  if (mfe>0){
	struc = space(strlen(seq)+1);
	E = fold(seq, struc);
	free(struc);
  } else {
	E=mfe;
  }
  
  avg = 0.0;
  stdv = 0.0;
  
  predict_values(seq, &avg, &stdv, type, avoid_shuffle, warning_string);

  /* Just as backup strategy if something goes totally wrong, 
     we evaluate the sequence once again by shuffling */
  if (avg > -1 || stdv < 0.1) {
    if (*type == 2) *type = 3;
    if (*type == 0) *type = 1;
    predict_values(seq, &avg, &stdv, type, avoid_shuffle, warning_string);
  }

  if ( stdv == 0.0) {
    return 0.0;
  }

  /*printf("%f,%f\n",avg,stdv);*/

  return ((E-avg)/stdv);

}
  
