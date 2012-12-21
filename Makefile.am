AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@

ACLOCAL_AMFLAGS="-Im4"

my_include_paths = -I$(top_srcdir)

AM_CXXFLAGS = $(CXXFLAGS) -I$(my_include_paths) $(BOOST_CPPFLAGS)
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
AM_FFLAGS = $(my_include_paths) -x f77-cpp-input

my_ld_flags = $(LDFLAGS) -lgsl -lgslcblas $(BOOST_PROGRAM_OPTIONS_LIBS) @FLIBS@
my_ld_paths = -L$(top_srcdir)  

# make reliefseq
if STATIC
bin_PROGRAMS = reliefseq reliefseq_static
else
bin_PROGRAMS = reliefseq
endif

reliefseq_SOURCES = ReliefSeqCLI.cpp Insilico.cpp DistanceMetrics.cpp \
Statistics.cpp Dataset.cpp ArffDataset.cpp \
PlinkDataset.cpp  PlinkBinaryDataset.cpp PlinkRawDataset.cpp DgeData.cpp \
BirdseedData.cpp DatasetInstance.cpp AttributeRanker.cpp ChiSquared.cpp \
ReliefF.cpp RReliefF.cpp SNReliefF.cpp ReliefFSeq.cpp ReliefSeqController.cpp

# libtool libraries
reliefseq_LDADD = $(my_ld_flags) $(my_ld_paths)

reliefseq_static_SOURCES=ReliefSeqCLI.cpp Insilico.cpp DistanceMetrics.cpp \
Statistics.cpp Dataset.cpp ArffDataset.cpp \
PlinkDataset.cpp  PlinkBinaryDataset.cpp PlinkRawDataset.cpp DgeData.cpp \
BirdseedData.cpp DatasetInstance.cpp AttributeRanker.cpp ChiSquared.cpp \
ReliefF.cpp RReliefF.cpp SNReliefF.cpp ReliefFSeq.cpp ReliefSeqController.cpp

# libtool libraries
reliefseq_static_LDFLAGS = -all-static -static-libgcc
reliefseq_static_LDADD = $(my_ld_flags) $(my_ld_paths)

# additional include paths to support the library headers
INCLUDES = $(myincludes)
