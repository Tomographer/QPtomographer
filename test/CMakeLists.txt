
#
# The tests can only be run with CMake for the moment.
#

# === Add the tests.

add_definitions(-DBOOST_TEST_DYN_LINK)

# test_utils:
add_executable(test_utils test_utils.cxx)
#set_property(TARGET test_utils PROPERTY CXX_STANDARD 11) # Enable C++11 -- flags already set by TomographerPy
target_include_directories(test_utils PUBLIC ${CMAKE_SOURCE_DIR}/cxx) # include our headers to test !
target_link_libraries(test_utils PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(test_utils PRIVATE TomographerPy::TomographerPy)
add_test(NAME test_utils
    COMMAND $<TARGET_FILE:test_utils>
    )

# test_diamond_norm:
add_executable(test_diamond_norm test_diamond_norm.cxx)
#set_property(TARGET test_diamond_norm PROPERTY CXX_STANDARD 11) # Enable C++11 -- flags already set by TomographerPy
target_link_libraries(test_diamond_norm PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_include_directories(test_diamond_norm PUBLIC ${CMAKE_SOURCE_DIR}/cxx)
target_include_directories(test_diamond_norm PUBLIC ${SDPA_INCLUDE_DIRS})
target_link_libraries(test_diamond_norm PRIVATE TomographerPy::TomographerPy)
target_link_libraries(test_diamond_norm PRIVATE ${SCS_DIR_LIBRARIES})
add_test(NAME test_diamond_norm
    COMMAND $<TARGET_FILE:test_diamond_norm>
    )

# test_channelspace:
add_executable(test_channelspace test_channelspace.cxx)
#set_property(TARGET test_channelspace PROPERTY CXX_STANDARD 11) # Enable C++11 -- flags already set by TomographerPy
target_include_directories(test_channelspace PUBLIC ${CMAKE_SOURCE_DIR}/cxx)
target_link_libraries(test_channelspace PRIVATE TomographerPy::TomographerPy)
target_link_libraries(test_channelspace PRIVATE ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_channelspace
    COMMAND $<TARGET_FILE:test_channelspace>
    )


#
# test that we can load and execute the python modules
#

add_test(NAME test_pyload_dnormtomo_bistates
  COMMAND ${PYTHON_EXECUTABLE} -c "import dnormtomo; print('dnormtomo: {}'.format(dnormtomo.__version__)); import sys; import numpy as np; import tomographer; import dnormtomo.bistates; dnormtomo.bistates.run(dimX=2,dimY=2,Emn=[np.diag([1,0,0,0])],Nm=[1000],hist_params=tomographer.HistogramParams(0,1,10),mhrw_params=tomographer.MHRWParams(0.001,1000,1024,4096),progress_fn=lambda x:sys.stdout.write(x.getHumanReport()+'\\n'))"
  )
set_tests_properties(test_pyload_dnormtomo_bistates PROPERTIES  ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR})

add_test(NAME test_pyload_dnormtomo_channelspace
  COMMAND ${PYTHON_EXECUTABLE} -c "import dnormtomo; print('dnormtomo: {}'.format(dnormtomo.__version__)); import sys; import numpy as np; import tomographer; import dnormtomo.channelspace; dnormtomo.channelspace.run(dimX=2,dimY=2,Emn=[np.diag([1,0,0,0])],Nm=[1000],hist_params=tomographer.HistogramParams(0,1,10),mhrw_params=tomographer.MHRWParams(0.001,1000,1024,4096),progress_fn=lambda x:sys.stdout.write(x.getHumanReport()+'\\n'))"
  )
set_tests_properties(test_pyload_dnormtomo_channelspace PROPERTIES  ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR})

