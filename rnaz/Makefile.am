AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I../librna -I../libsvm-@LIBSVM_VERSION@ -I../models

LDADD = -L../librna -lRNA -lm
 
bin_PROGRAMS = RNAz

RNAz_SOURCES=RNAz.c rnaz_utils.c svm_helper.c rnaz_utils.h svm_helper.h zscore.c zscore.h cmdline.c cmdline.h strand.c strand.h  ../models/mfe_avg.inc  ../models/mfe_stdv.inc  ../models/strand.inc ../models/decision.inc 

nodist_RNAz_SOURCES=../libsvm-@LIBSVM_VERSION@/svm.cpp 

RNAz_LINK = $(CXX) -o $@
 
