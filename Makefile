HD_COMPONENT_NAME	= Xspec

HD_COMPONENT_VERS	= 

#library definition

HD_LIBRARY_ROOT		= mcvspec

HD_LIB_STYLE		= shared

#source code

HD_LIBRARY_SRC_f	= 


HD_LIBRARY_SRC_f90	= mcvspec.f90


HD_LIBRARY_SRC_f03	= 

HD_LIBRARY_SRC_c	= 

HD_LIBRARY_SRC_C	= 

HD_LIBRARY_SRC_cxx	= lpack_mcvspec.cxx \
			  mcvspecFunctionMap.cxx

HD_LIBRARY_SRC_cpp	= 

HD_LIBRARY_SRC_cc	= 


HD_INSTALL_LIBRARIES	= ${HD_LIBRARY_ROOT}

HD_CXXFLAGS		= ${HD_STD_CXXFLAGS} \
			  -I${HEADAS}/include -DINITPACKAGE

HD_CFLAGS		= ${HD_STD_CFLAGS} \
			  -I${HEADAS}/include -DINITPACKAGE

HD_FFLAGS		= ${HD_STD_FFLAGS} \
			  -I${HEADAS}/include -DINITPACKAGE

#lib file name
PACKAGE		= lib${HD_LIBRARY_ROOT}${SHLIB_SUFFIX}

HD_CLEAN		= lpack_${HD_LIBRARY_ROOT}.cxx \
			  ${HD_LIBRARY_ROOT}FunctionMap.cxx \
			  ${HD_LIBRARY_ROOT}FunctionMap.h \
			  ${PACKAGE} Makefile pkgIndex.tcl *.bck

HD_SHLIB_LIBS           = ${HD_LFLAGS} -l${CCFITS} -l${CFITSIO} -lXS \
                          -lXSUtil -lXSFunctions -lXSModel -l${TCL} \
                          -l${HEAUTILS} ${HD_STD_LIBS} ${SYSLIBS} ${F77LIBS4C} \

include ${HD_STD_MAKEFILE}

