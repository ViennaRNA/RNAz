bin_PROGRAMS = RNAz

RNAz_SOURCES=RNAz.c svm_helper.c svm_helper.h zscore.c zscore.h

LDADD = -L../librna -lRNA -L../libsvm -lsvm -lm
 
RNAz_LINK = $(CXX) -o $@
 
INCLUDES = -I../librna -I../libsvm

