TARGET := reGeneration_nGd_gamma

#SRCS = $(TARGET).cc
SRCDIR = src
LIBDIR = lib
SRCS = $(wildcard $(SRCDIR)/*.cc) $(wildcard $(LIBDIR)/*.cc)
OBJDIR = obj
#OBJS = $(TARGET).o
#OBJS = $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cc=.o)))
OBJS = $(addprefix $(OBJDIR)/, $(patsubst %.cc, %.o, ${SRCS}))

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
JADEDIR = /group/mlf/nu/work/software/jade
JADELIBS = -L$(JADEDIR)/lib
JADEINC = -I$(JADEDIR)/include
RATINC = -I$(RATROOT)/include
RATLIBS = -L$(RATROOT)/lib -lRATEvent -lrat

ANNINC = include

CXXFLAGS = $(ROOTCFLAGS) -Wall $(JADEINC) $(RATINC) -I$(ANNINC)
CXXLIBS = $(ROOTLIBS) $(ROOTGLIBS) $(JADELIBS) $(RATLIBS) -lMinuit 
CXX = g++

all: $(TARGET)

$(OBJDIR)/%.o : ./%.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<
#.cc.o:
#	$(CXX) $(CXXFLAGS) $(JADELIBS) -c $<

$(TARGET): $(OBJS)
	$(CXX) $(CXXLIBS) -o $@ $^

clean:
	rm -f $(TARGET)
