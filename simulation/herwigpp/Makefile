# -*- Makefile -*- (for emacs)
# src/Makefile-UserModules.  Generated from Makefile-UserModules.in by configure.
#
# Herwig++ 2.7.1
#
# This Makefile is intended for compiling Herwig++ plugins
# You can find example plugins in the Contrib directory
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#

#
# C++ flags
# 
CXX = g++

CXXFLAGS = -g -W -Wall -O3

# add additional package header directories here
CPPFLAGS = -I/usr/local/include -I/usr/local/include 

# add additional libraries here
LDFLAGS = 

ALLCCFILES=$(shell echo *.cc)
all : $(ALLCCFILES:.cc=.so)

%.so : %.cc %.h 
	$(CXX) -shared -fPIC $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< -o $@  

clean:
	rm -f $(ALLCCFILES:.cc=.so)
