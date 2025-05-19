# $Id: GNUmakefile,v 1.2 2000/10/19 12:22:10 stanaka Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := exe_ANNRI_BGO
G4TARGET := $(name)
G4EXLIB := true
G4DEBUG := true
G4PROFILE := true

G4WORKDIR := .

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
#EXTRALIBS = -pg
EXTRALIBS =

ifdef ROOTSYS
  ifeq ($(G4SYSTEM), Linux-g++)
    #  for Linux
    CPPFLAGS += -g $(EXTRALIBS) -O2 -I$(shell $(ROOTSYS)/bin/root-config --incdir)
    #LDLIBS   += -g -Wl,-rpath,$(shell $(ROOTSYS)/bin/root-config --libdir) 
    LDLIBS   += -Wl,-rpath,$(shell $(ROOTSYS)/bin/root-config --libdir) \
                $(shell $(ROOTSYS)/bin/root-config --libs)
  else
    CPPFLAGS += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
    LDLIBS   += $(shell $(ROOTSYS)/bin/root-config --libs)
  endif
endif

