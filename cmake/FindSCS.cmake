
# Find SCS library
#
# This script defines the following variables:
#
#  SCS_FOUND - SCS was found on this system
#  SCS_INCLUDE_DIRS - include directories (for -I flags)
#  SCS_DIR_LIBRARIES - libraries needed to link to SCS
#  SCS_INDIR_LIBRARIES - libraries needed to link to SCS
#  SCS_DEFINITIONS - Compiler definitions needed for SCS
#
# You may specify hints as:
#
#  SCS_ROOT - directory of downloaded & compiled scs library
#  SCS_DONT_PREFER_STATIC - don't prefer static libraries over the shared libraries
#

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if (SCS_DONT_PREFER_STATIC)
  set(_scs_staticlib_name_dir )
  set(_scs_staticlib_name_indir )
else()
  set(_scs_staticlib_name_dir   libscsdir.a)
  set(_scs_staticlib_name_indir libscsindir.a)
endif()

find_path(SCS_INCLUDE_DIR
  scs.h
  HINTS ${SCS_ROOT}/include
  )
find_path(SCSLINSYS_INCLUDE_DIR
  linsys/amatrix.h
  HINTS ${SCS_ROOT}
  )
find_library(SCS_DIR_LIBRARY
  NAMES ${_scs_staticlib_name_dir} scsdir
  HINTS ${SCS_ROOT}/out
  )
find_library(SCS_INDIR_LIBRARY
  NAMES ${_scs_staticlib_name_indir} scsindir
  HINTS ${SCS_ROOT}/out
  )

set(SCS_DIR_LIBRARIES ${SCS_DIR_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set(SCS_INDIR_LIBRARIES ${SCS_INDIR_LIBRARY}  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set(SCS_INCLUDE_DIRS ${SCS_INCLUDE_DIR} ${SCSLINSYS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SCS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SCS  DEFAULT_MSG
                                  SCS_DIR_LIBRARY SCS_INDIR_LIBRARY SCS_INCLUDE_DIR SCSLINSYS_INCLUDE_DIR)

mark_as_advanced(SCS_INCLUDE_DIR SCSLINSYS_INCLUDE_DIR SCS_DIR_LIBRARY SCS_INDIR_LIBRARY )
