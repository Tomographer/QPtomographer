
#
# will find also Python libs & interpreter
#
#

find_package(PythonLibsNew)

# No, use Python's pybind11.
#
#find_package(pybind11 2.1 REQUIRED)


# write a small python script that will find everything & generate cmake code that we can import
set(_tomographerpy_cmake_setcode "
import sys; import os; import re; import tomographer.include; import pybind11;
def q(x):
    return '\"' + re.sub(r'(\"|\\\\)', lambda m: r'\\\\'+m.group(), x) + '\"'
if not sys.argv[1]:
    raise RuntimeError('Invalid {fn}!'.format(fn=sys.argv[1]))
f = open(sys.argv[1], 'w')
def print_set(var, x):
    f.write('set('+var+'  '+q(x)+')\\n')
f.write('macro(SetupTomographerPyVars)\\n')
print_set('TomographerPy_INCLUDE_DIR_TOMOGRAPHER', tomographer.include.get_include(keys=True)['tomographer'])
print_set('TomographerPy_INCLUDE_DIR_TOMOGRAPHERPY', tomographer.include.get_include(keys=True)['tomographerpy'])
print_set('TomographerPy_INCLUDE_DIR_DEPS',
          ';'.join([v for k,v in tomographer.include.get_include(keys=True).items() if k != 'tomographerpy' and k != 'tomographer']))
print_set('TomographerPy_INCLUDE_DIR_DEP_EIGEN', tomographer.include.get_include(keys=True)['eigen'])
print_set('TomographerPy_INCLUDE_DIR_DEP_BOOST', tomographer.include.get_include(keys=True)['boost'])
print_set('TomographerPy_INCLUDE_DIR_PYBIND11', pybind11.get_include())
print_set('TomographerPy_CXX_FLAGS', ';'.join(tomographer.version.compile_info['cflags']))
print_set('TomographerPy_VERSION_STR', tomographer.__version__)
print_set('TomographerPy_VERSION', '.'.join([str(x) for x in tomographer.version.version_info]))
f.write('endmacro(SetupTomographerPyVars)\\n')
f.close()
")
file(WRITE "${CMAKE_BINARY_DIR}/tmp_find_tomographer_py_mk_cmds.py" "${_tomographerpy_cmake_setcode}")

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_BINARY_DIR}/tmp_find_tomographer_py_mk_cmds.py" "${CMAKE_BINARY_DIR}/tmp_find_tomographer_py_cmds.cmake"
  RESULT_VARIABLE _tomographerpy_result
  #OUTPUT_VARIABLE _tomographerpy_setcode
  #OUTPUT_STRIP_TRAILING_WHITESPACE
  )
if(_tomographerpy_result)
  message(FATAL_ERROR "Couldn't find TomographerPy: Error: ${_tomographerpy_result}")
endif()

include("${CMAKE_BINARY_DIR}/tmp_find_tomographer_py_cmds.cmake")
SetupTomographerPyVars()


add_library(TomographerPy::TomographerRaw INTERFACE IMPORTED)
set_target_properties(TomographerPy::TomographerRaw  PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TomographerPy_INCLUDE_DIR_TOMOGRAPHER}"
  )
add_library(TomographerPy::Deps::Eigen INTERFACE IMPORTED)
set_target_properties(TomographerPy::Deps::Eigen  PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TomographerPy_INCLUDE_DIR_DEP_EIGEN}"
  )
add_library(TomographerPy::Deps::BoostHeaders INTERFACE IMPORTED)
set_target_properties(TomographerPy::Deps::BoostHeaders  PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TomographerPy_INCLUDE_DIR_DEP_BOOST}"
  )


add_library(TomographerPy::Tomographer INTERFACE IMPORTED)
set_target_properties(TomographerPy::Tomographer  PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TomographerPy_INCLUDE_DIR_DEPS}"
  )
add_dependencies(TomographerPy::Tomographer TomographerPy::TomographerRaw)

add_library(TomographerPy::TomographerPy INTERFACE IMPORTED)
set_target_properties(TomographerPy::TomographerPy  PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TomographerPy_INCLUDE_DIR_TOMOGRAPHERPY};${PYTHON_INCLUDE_DIRS}"
  INTERFACE_COMPILE_OPTIONS "${TomographerPy_CXX_FLAGS}"
  )
add_dependencies(TomographerPy::TomographerPy TomographerPy::Tomographer)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TomographerPy
  FOUND_VAR TomographerPy_FOUND
  REQUIRED_VARS
    TomographerPy_INCLUDE_DIR_TOMOGRAPHER
    TomographerPy_INCLUDE_DIR_TOMOGRAPHERPY
    TomographerPy_INCLUDE_DIR_DEP_EIGEN
    TomographerPy_INCLUDE_DIR_DEP_BOOST
    TomographerPy_INCLUDE_DIR_DEPS
    TomographerPy_CXX_FLAGS
    TomographerPy_VERSION_STR
    TomographerPy_VERSION
  VERSION_VAR TomographerPy_VERSION
  )


macro(SetupTomographerPyModule target_name)
  target_link_libraries(${target_name} PRIVATE TomographerPy::TomographerPy)
  set_target_properties(${target_name} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
  set_target_properties(${target_name} PROPERTIES SUFFIX "${PYTHON_MODULE_EXTENSION}")
  target_link_libraries(${target_name} PRIVATE "-undefined dynamic_lookup")
  # will fall back to an earlier standard if C++14 isn't supported
  set_target_properties(${target_name} PROPERTIES CXX_STANDARD 14)
endmacro()

