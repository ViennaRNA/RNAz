# Inits project from configure.ac
aclocal
touch stamp-h
autoheader
autoconf
automake -a
