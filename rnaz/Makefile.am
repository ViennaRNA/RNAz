bin_PROGRAMS = RNAz

RNAz_SOURCES=RNAz.c rnaz_utils.c svm_helper.c rnaz_utils.h svm_helper.h zscore.c zscore.h cmdline.c cmdline.h strand.c strand.h  ../models/mfe_avg.inc  ../models/mfe_stdv.inc  ../models/strand.inc ../models/decision.inc 

nodist_RNAz_SOURCES=../libsvm-2.8/svm.cpp 

LDADD = -L../librna -lRNA -lm
 
RNAz_LINK = $(CXX) -o $@
 
INCLUDES = -I../librna -I../libsvm-2.8 -I../models

