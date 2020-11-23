IDIR=include

CXX=clang++
CXXFLAGS=-I$(IDIR) -std=c++2a -Wall -Wextra -Wpedantic

SDIR=src
ODIR=$(SDIR)/obj

_SRCS=main.cpp kernels.cpp eos.cpp vtkout.cpp sph.cpp
SRCS=$(patsubst %,$(SDIR)/%,$(_SRCS))

_OBJS=$(_SRCS:.cpp=.o)
OBJS=$(patsubst %,$(ODIR)/%,$(_OBJS))

_DEPS=kernels.h eos.h particle.h vtkout.h sph.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

EXE=sph-lite.exe
DBGEXE=sph-lite.exe.debug

all: release

release: CXXFLAGS += -O3
release: clean $(EXE)

parallel: CXXFLAGS += -fopenmp
parallel: clean release

debug: CXXFLAGS += -g
debug: clean $(DBGEXE)

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(EXE) $(DBGEXE): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(DBGEXE) $(EXE)
