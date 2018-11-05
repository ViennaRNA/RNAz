AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = \
    -I$(top_srcdir)/librna \
    -I$(top_srcdir)/libsvm-@LIBSVM_VERSION@ \
    -I$(top_srcdir)/models

LDADD = \
    -L$(top_builddir)/librna \
    -lRNA \
    -lm
 
bin_PROGRAMS = RNAz

noinst_HEADERS = \
    rnaz_utils.h \
    svm_helper.h \
    zscore.h \
    cmdline.h \
    strand.h

SVM_MODEL_INC = \
    $(top_srcdir)/models/mfe_avg.inc \
    $(top_srcdir)/models/mfe_stdv.inc \
    $(top_srcdir)/models/strand.inc \
    $(top_srcdir)/models/decision.inc


RNAz_SOURCES = \
    RNAz.c \
    rnaz_utils.c \
    svm_helper.c \
    zscore.c \
    cmdline.c \
    strand.c \
    ../libsvm-@LIBSVM_VERSION@/svm.cpp


RNAz_LINK = $(CXX) -o $@


EXTRA_DIST = \
    ${SVM_MODEL_INC}
