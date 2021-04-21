BINDIR   := bin
CXXFLAGS := -g -O2 -Wall -Wextra -std=c++17 -pedantic $(CXXFLAGS)
LIBS     :=

# Targets
EXE      := $(BINDIR)/bditau
EXESRC   := $(patsubst $(BINDIR)/%,%.cc,$(EXE))
EXEOBJ   := $(EXESRC:.cc=.o)

# ROOT (https://root.cern/)
CXXFLAGS += $(shell root-config --cflags)
LIBS     += $(shell root-config --libs)

# YAM2 (https://github.com/cbpark/YAM2)
YAM2     ?= /usr/local
CXXFLAGS += -I$(YAM2)/include
LIBS     += -L$(YAM2)/lib -lYAM2 -Wl,-rpath $(YAM2)/lib

# NLopt (https://nlopt.readthedocs.io/)
NLOPT    ?= /usr
LIBS     += -L$(NLOPT)/lib -lnlopt -Wl,-rpath $(NLOPT)/lib

.PHONY: all build clean

all: $(EXE)

$(EXE): $(EXEOBJ) build
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

build:
	mkdir -p $(BINDIR)

clean:
	rm -f $(EXE) $(EXEOBJ)
	rmdir $(BINDIR)
