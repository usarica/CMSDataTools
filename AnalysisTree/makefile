COMPILEPATH          = $(PWD)/
BASEINCLUDE          = $(COMPILEPATH)/../../

INCLUDEDIR           = $(COMPILEPATH)interface/
SRCDIR               = $(COMPILEPATH)src/
BINDIR               = $(COMPILEPATH)bin/
LIBDIR               = $(COMPILEPATH)lib/
TESTDIR              = $(COMPILEPATH)test/
PYTHONDIR            = $(COMPILEPATH)python/
RUNDIR               = $(COMPILEPATH)
LIB                  = libCMSDataToolsAnalysisTree.so

MELADIR              = $(BASEINCLUDE)JHUGenMELA/MELA
MELALIBDIR           = $(MELADIR)/data/${SCRAM_ARCH}

GMECOMPDIR           = $(BASEINCLUDE)MelaAnalytics/GenericMEComputer
GMECOMPLIBDIR        = $(GMECOMPDIR)/lib
GMECOMPLIBNAME       = libMelaAnalyticsGenericMEComputer.so

MELAEVTCONTDIR       = $(BASEINCLUDE)MelaAnalytics/EventContainer
MELAEVTCONTLIBDIR    = $(MELAEVTCONT)/lib
MELAEVTCONTLIBNAME   = libMelaAnalyticsEventContainer.so

MELARECASTERDIR      = $(BASEINCLUDE)MelaAnalytics/CandidateLOCaster
MELARECASTERLIBDIR   = $(MELARECASTERDIR)/lib
MELARECASTERLIBNAME  = libMelaAnalyticsCandidateLOCaster.so


MELACXXFLAGS =  -I$(MELADIR)/interface/ -L$(MELALIBDIR)
MELALIBS =  -lJHUGenMELAMELA

MELAANALYTICSLIBS = -lMelaAnalyticsGenericMEComputer -lMelaAnalyticsEventContainer -lMelaAnalyticsCandidateLOCaster
MELAANALYTICSCXXFLAGS = -I$(GMECOMPDIR)/interface/ -L$(GMECOMPLIBDIR) -I$(MELAEVTCONTDIR)/interface/ -L$(MELAEVTCONTLIBDIR) -I$(MELARECASTERDIR)/interface/ -L$(MELARECASTERLIBDIR)


EXTCXXFLAGS   = $(MELACXXFLAGS)
EXTLIBS       = $(MELALIBS)

ROOTCFLAGS    = $(shell root-config --cflags) -Lrootlib
ROOTLIBS     = $(shell root-config --libs) -Lrootlib

ARCH         := $(shell root-config --arch)

CXX           = g++
CXXFLAGS      = -fPIC -g -O2 $(ROOTCFLAGS) -D_COMPILE_STANDALONE_ -I$(ROOFITSYS)/include/ -I$(INCLUDEDIR) $(EXTCXXFLAGS)
LINKERFLAGS   = -Wl,-rpath=$(LIBDIR),-soname,$(LIB)

NLIBS         = $(ROOTLIBS)
# Hack here, because RooFit is removed from ROOT:
NLIBS        += -L $(ROOFITSYS)/lib/ -lMinuit -lRooFitCore -lRooFit #-lboost_regex
# Libraries for common user packages
NLIBS        += $(EXTLIBS)
# Filter out libNew because it leads to floating-point exceptions in versions of ROOT prior to 6.08.02
# See the thread https://root-forum.cern.ch/t/linking-new-library-leads-to-a-floating-point-exception-at-startup/22404
LIBS          = $(filter-out -lNew, $(NLIBS))


SOURCES = $(wildcard $(SRCDIR)*.cc)
OBJECTS = $(SOURCES:.cc=.o)

.PHONY: all help compile clean

all: $(LIB) python/__init__.py MergeLHEFiles CheckFileIntegrity SplitROOTFiles ExecuteCompiledCommand


# General rule for making object files
%.d:	$(SRCDIR)%.cc $(BINDIR)%.cc $(SRCDIR)%.cxx $(BINDIR)%.cxx
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CXX) -MM -MT $@ -MT ${@:.d=.o} $(CXXFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	$(SRCDIR)%.cc $(SRCDIR)%.cxx $(BINDIR)%.cc $(BINDIR)%.cxx
	$(QUIET) echo "Compiling $<"; \
	$(CXX) $(CXXFLAGS) $< -c -o $@ $(LIBS)


$(LIB):	$(OBJECTS)
	$(QUIET) echo "Linking $(LIB)"; \
	$(CXX) $(LINKERFLAGS) -shared $(OBJECTS) -o $@; \
        mkdir -p $(LIBDIR); \
	mv $(LIB) $(LIBDIR)/$(LIB)

python/__init__.py:
	touch $(PYTHONDIR)__init__.py

MergeLHEFiles: $(LIB) $(BINDIR)MergeLHEFiles.cc
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)MergeLHEFiles $(BINDIR)MergeLHEFiles.cc $(LIBS) -L$(LIBDIR) -lCMSDataToolsAnalysisTree

CheckFileIntegrity: $(LIB) $(BINDIR)CheckFileIntegrity.cc
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)CheckFileIntegrity $(OBJECTS) $(BINDIR)CheckFileIntegrity.cc $(LIBS) -L$(LIBDIR) -lCMSDataToolsAnalysisTree

SplitROOTFiles: $(LIB) $(BINDIR)SplitROOTFiles.cc
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)SplitROOTFiles $(OBJECTS) $(BINDIR)SplitROOTFiles.cc $(LIBS) -L$(LIBDIR) -lCMSDataToolsAnalysisTree

ExecuteCompiledCommand: $(LIB) $(BINDIR)ExecuteCompiledCommand.cc
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)ExecuteCompiledCommand $(OBJECTS) $(BINDIR)ExecuteCompiledCommand.cc $(LIBS) -L$(LIBDIR) -lCMSDataToolsAnalysisTree

clean:
	rm -rf $(LIBDIR)
	rm -f $(SRCDIR)*.o
	rm -f $(SRCDIR)*.so
	rm -f $(SRCDIR)*.d
	rm -f $(BINDIR)*.o
	rm -f $(BINDIR)*.so
	rm -f $(BINDIR)*.d
	rm -f $(PYTHONDIR)__init__.py
	rm -f $(PYTHONDIR)*.pyc
	rm -rf $(RUNDIR)Pdfdata
	rm -f $(RUNDIR)input.DAT
	rm -f $(RUNDIR)process.DAT
	rm -f $(RUNDIR)ffwarn.dat
	rm -f $(RUNDIR)br.sm*
	rm -f $(RUNDIR)*.cc
	rm -f $(RUNDIR)*.o
	rm -f $(RUNDIR)*.so
	rm -f $(RUNDIR)*.d
	rm -f $(RUNDIR)*.pcm
	rm -f $(RUNDIR)*.pyc
	rm -f $(RUNDIR)MergeLHEFiles
	rm -f $(RUNDIR)CheckFileIntegrity
	rm -f $(RUNDIR)SplitROOTFiles
	rm -f $(RUNDIR)ExecuteCompiledCommand
	rm -rf $(TESTDIR)Pdfdata
	rm -f $(TESTDIR)input.DAT
	rm -f $(TESTDIR)process.DAT
	rm -f $(TESTDIR)ffwarn.dat
	rm -f $(TESTDIR)br.sm*
	rm -f $(TESTDIR)*.o
	rm -f $(TESTDIR)*.so
	rm -f $(TESTDIR)*.d
	rm -f $(TESTDIR)*.pcm
	rm -f $(TESTDIR)*.pyc

