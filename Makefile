# This is probably painful for anyone who can actually write makefiles
CXXFLAGS  = -std=c++17
# CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -fPIC
CXXFLAGS += -pipe
CXXFLAGS += -D_VANILLA_ROOT_
CXXFLAGS += -g
# CXXFLAGS += -no-pie
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

qa_object_list = qa.o setup.o
post_object_list = qa_histograms.o setup.o draw_histogram.o histogram_package.o histogram_data.o
analysis_object_list = calculate_v2.o setup.o histogram_package.o histogram_data.o

qa_objects = $(qa_object_list:%.o=build/%.o)
post_objects = $(post_object_list:%.o=build/%.o)
analysis_objects = $(analysis_object_list:%.o=build/%.o)


all: qa post_qa analysis

build/%.o: %.cxx 
	$(CXX) -c $(CXXFLAGS) $(INCFLAGS) $< -o $@

qa: $(qa_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBPATH) $(qa_objects) $(LIBS) -o qa

post_qa: $(post_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBPATH) $(post_objects) $(LIBS) -o post_qa

analysis: $(analysis_objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBPATH) $(analysis_objects) $(LIBS) -o analysis


clean:
	rm -f build/*
	rm -f qa
	rm -f post_qa
	rm -f analysis

