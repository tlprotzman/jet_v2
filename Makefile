CXXFLAGS  = -std=c++17
# CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -fPIC
CXXFLAGS += -pipe
CXXFLAGS += -D_VANILLA_ROOT_
CXXFLAGS += -g

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
LIBS += -lStPicoEvent
LIBS += -lTStarJetPico
LIBS += -ljetreader
LIBS += -lStEpdUtil

all:
	$(CXX) -c $(CXXFLAGS) $(INCFLAGS) qa.cxx -o qa.o
	$(CXX) $(LDFLAGS) $(LIBPATH) qa.o $(LIBS) -o qa.out