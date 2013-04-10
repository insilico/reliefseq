AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@

ACLOCAL_AMFLAGS="-Im4"

AM_CXXFLAGS = -I$(top_srcdir) $(BOOST_CPPFLAGS)
# AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
AM_FFLAGS = $(my_include_paths)

my_ld_flags = $(LDFLAGS) -lgsl -lgslcblas $(BOOST_PROGRAM_OPTIONS_LIBS)
my_ld_paths = -L$(top_srcdir)  

bin_PROGRAMS = reliefseq

reliefseq_SOURCES = ReliefSeqCLI.cpp Insilico.cpp DistanceMetrics.cpp \
Statistics.cpp Dataset.cpp ArffDataset.cpp \
PlinkDataset.cpp  PlinkBinaryDataset.cpp PlinkRawDataset.cpp DgeData.cpp \
BirdseedData.cpp DatasetInstance.cpp AttributeRanker.cpp ChiSquared.cpp \
ReliefF.cpp RReliefF.cpp SNReliefF.cpp ReliefFSeq.cpp ReliefSeqController.cpp \
config.h GSLRandomBase.h GSLRandomFlat.h Insilico.h DistanceMetrics.h \
Statistics.h Dataset.h ArffDataset.h StringUtils.h BestN.h \
PlinkDataset.h  PlinkBinaryDataset.h PlinkRawDataset.h DgeData.h \
BirdseedData.h DatasetInstance.h AttributeRanker.h ChiSquared.h \
ReliefF.h RReliefF.h SNReliefF.h ReliefFSeq.h ReliefSeqController.h

# libtool libraries
reliefseq_LDFLAGS = -fopenmp
reliefseq_LDADD = $(my_ld_flags) $(my_ld_paths)

# additional include paths to support the library headers
INCLUDES = $(myincludes)
