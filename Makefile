CXXFLAGS  = -std=c++17
# CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -fPIC
CXXFLAGS += -pipe
CXXFLAGS += -D_VANILLA_ROOT_
CXXFLAGS += -g
CXXFLAGS += -O3

INCFLAGS  = -I$(ROOTSYS)/include
INCFLAGS += -I$(FASTJETDIR)/include
INCFLAGS += -I$(STARPICOPATH)
INCFLAGS += -I$(JETREADER)/include

LIBPATH  = $(shell root-config --libs)
LIBPATH += -L$(FASTJETDIR)/lib
LIBPATH += -L$(JETREADER)/lib
LIBPATH += -L$(STARPICOPATH)
LIBPATH += -LStEpdUtil

LIBS  = -lfastjet
LIBS += -lfastjettools
LIBS += -lStPicoEvent
LIBS += -lTStarJetPico
LIBS += -ljetreader
LIBS += -lStEpdUtil
LIBS += -lNetx

objects = qa.o setup.o

%.o: %.cxx
	$(CXX) -c $(CXXFLAGS) $(INCFLAGS) $< -o $@

qa.out: $(objects)
	$(CXX) $(LDFLAGS) $(LIBPATH) $(objects) $(LIBS) -o qa.out

all: qa.out
