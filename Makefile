IDIR=include

CXX=g++
CXXFLAGS=-I$(IDIR) -O3

SDIR=src
ODIR=$(SDIR)/obj

_SRCS=sph-lite.cpp kernels.cpp eos.cpp
SRCS=$(patsubst %,$(SDIR)/%,$(_SRCS))

_OBJS=$(_SRCS:.cpp=.o)
OBJS=$(patsubst %,$(ODIR)/%,$(_OBJS))

_DEPS=kernels.h eos.h particle.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

EXE=sph-lite.exe
DBGEXE=sph-lite.exe.debug

all: clean $(EXE)

debug: CXXFLAGS += -g
debug: $(DBGEXE)

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(EXE) $(DBGEXE): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(DBGEXE) $(EXE)
