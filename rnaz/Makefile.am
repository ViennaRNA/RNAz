bin_PROGRAMS = RNAz

RNAz_SOURCES=RNAz.c svm_helper.c svm_helper.h zscore.c zscore.h 

nodist_RNAz_SOURCES=../libsvm-2.6/svm.cpp

LDADD = -L../librna -lRNA -lm
 
RNAz_LINK = $(CXX) -o $@
 
INCLUDES = -I../librna -I../libsvm-2.6

