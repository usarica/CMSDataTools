ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else
ifndef ROOFITSYS
all:
	@echo "ROOFITSYS is not set. Please set ROOT environment properly"; echo
else


COMPILEPATH          = $(PWD)/
BASEINCLUDE          = $(COMPILEPATH)../../

INCLUDEDIR           = $(COMPILEPATH)interface/
SRCDIR               = $(COMPILEPATH)src/
BINDIR               = $(COMPILEPATH)bin/
SCRIPTSDIR           = $(COMPILEPATH)scripts/
OBJDIR               = $(COMPILEPATH)obj/
LIBDIR               = $(COMPILEPATH)lib/
EXEDIR               = $(COMPILEPATH)executables/
TESTDIR              = $(COMPILEPATH)test/
PYTHONDIR            = $(COMPILEPATH)python/
RUNDIR               = $(COMPILEPATH)
LIB                  = libCMSDataToolsAnalysisTree.so
LIBRULE              = $(LIBDIR)$(LIB)

MELADIR              = $(BASEINCLUDE)JHUGenMELA/MELA/
MELALIBDIR           = $(MELADIR)data/${SCRAM_ARCH}


MELACXXFLAGS =  -I$(MELADIR)interface/ -L$(MELALIBDIR)
MELALIBS =  -lJHUGenMELAMELA


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


SOURCESCC = $(wildcard $(SRCDIR)*.cc)
SOURCESCXX = $(wildcard $(SRCDIR)*.cxx)
OBJECTSPRIM = $(SOURCESCC:.cc=.o) $(SOURCESCXX:.cxx=.o)
OBJECTS = $(subst $(SRCDIR),$(OBJDIR),$(OBJECTSPRIM))

BINSCC = $(wildcard $(BINDIR)*.cc)
BINSCXX = $(wildcard $(BINDIR)*.cxx)
EXESPRIM = $(BINSCC:.cc=) $(BINSCXX:.cxx=)
EXES = $(subst $(BINDIR),$(EXEDIR),$(EXESPRIM))


.PHONY: all help compile clean
.SILENT: alldirs scripts clean $(OBJECTS) $(LIBRULE) $(EXES)

all: $(OBJECTS) $(LIBRULE) python $(EXES) scripts


alldirs:
	mkdir -p $(OBJDIR); \
	mkdir -p $(LIBDIR); \
	mkdir -p $(EXEDIR)

scripts: | alldirs
	cp $(SCRIPTSDIR)* $(EXEDIR)

python:
	touch $(PYTHONDIR)__init__.py

$(LIBRULE):	$(OBJECTS) | alldirs
	echo "Linking $(LIB)"; \
	$(CXX) $(LINKERFLAGS) -shared $(OBJECTS) -o $@

$(OBJDIR)%.d:	$(SRCDIR)%.c* | alldirs
	echo "Checking dependencies for $<"; \
	$(CXX) -MM -MT $@ -MT ${@:.d=.o} $(CXXFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

$(OBJDIR)%.o: 	$(SRCDIR)%.c* | alldirs
	echo "Compiling $<"; \
	$(CXX) $(CXXFLAGS) $< -c -o $@ $(LIBS)


$(EXEDIR)%::	$(BINDIR)%.cc $(LIBRULE) | alldirs
	echo "Compiling $<"; \
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS) -L$(LIBDIR) -lCMSDataToolsAnalysisTree


clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(EXEDIR)
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



endif
endif
