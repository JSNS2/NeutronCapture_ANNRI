TARGET := reGeneration_nGd_gamma

#SRCS = $(TARGET).cc
INCDIR := include
SRCDIR := src
LIBDIR := lib
OUTDIR := obj
TARGET_2 = $(OUTDIR)/$(TARGET)
SRCS = $(wildcard $(SRCDIR)/*.cc) $(wildcard *.cc)
#SRCS = $(wildcard $(SRCDIR)/*.cc) $(wildcard $(LIBDIR)/*.cc)
#SRCS = $(wildcard $(SRCDIR)/*.cc)
OBJS = $(addprefix $(OUTDIR)/, $(patsubst %.cc, %.o, ${SRCS}))

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
RATINC = -I$(RATROOT)/include
RATLIBS = -L$(RATROOT)/lib -lRATEvent -lrat


CXXFLAGS = $(ROOTCFLAGS) -Wall $(RATINC) -I$(INCDIR)
CXXLIBS = $(ROOTLIBS) $(ROOTGLIBS) $(RATLIBS) -lMinuit 
CXX = g++

.PHONY: all clean
all: $(TARGET_2)

$(TARGET_2): $(OBJS)
	$(CXX) $(CXXLIBS) -o $@ $^

#$(OBJDIR)/%.o : ./%.cc
$(OUTDIR)/%.o : ./%.cc
	@if [ ! -e `dirname $@` ] ; then mkdir -p `dirname $@`; fi
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf $(OUTDIR)
