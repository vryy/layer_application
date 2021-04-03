set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosLayerApplication ****")

######### add include directory here
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
include_directories( ${CMAKE_SOURCE_DIR} )
#######################################################################

# zlib, lib64

set( libb64_src
    ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/libb64/cdecode.c
    ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/libb64/cencode.c
)
add_library(b64 SHARED ${libb64_src})
install(TARGETS b64 DESTINATION libs)

add_definitions(-D_USE_ZLIB)
add_definitions(-D_USE_LIBB64)

# hdf5
if(DEFINED $ENV{HDF5_ROOT})
    # SET(HDF5_DIR $ENV{HDF5_ROOT}/share/cmake/hdf5)
    # FIND_PACKAGE(HDF5 QUIET NAMES hdf5 hdf5_cpp)
    FIND_PACKAGE(HDF5 REQUIRED COMPONENTS CXX)
    set(HDF5_PROVIDED TRUE)
elseif(DEFINED HDF5_ROOT)
    # SET(HDF5_DIR ${HDF5_ROOT}/share/cmake/hdf5)
    FIND_PACKAGE(HDF5 REQUIRED COMPONENTS CXX)
    set(HDF5_PROVIDED TRUE)
else()
    set(HDF5_PROVIDED FALSE)
endif()

set(LAYER_APP_USE_HDF5 FALSE)
if(${HDF5_PROVIDED} MATCHES TRUE)
    message("HDF5_INCLUDE_DIR: " ${HDF5_INCLUDE_DIR})
    message("HDF5_LIBRARIES: " ${HDF5_LIBRARIES})
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
#    add_definitions( -DBOOST_FILESYSTEM_VERSION=2 )
    add_definitions(-DLAYER_APP_USE_HDF5)
    set(LAYER_APP_USE_HDF5 TRUE)
endif()

######### add compilation units here
set( KRATOS_LAYER_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/post_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/model_part_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/layer_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_io_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_parameter_list_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/layer_python_application.cpp
)
#######################################################################

## define library Kratos which defines the basic python interface
add_library(KratosLayerApplication SHARED ${KRATOS_LAYER_APPLICATION_SOURCES})
target_link_libraries(KratosLayerApplication KratosCore)
target_link_libraries(KratosLayerApplication zlibstatic)
target_link_libraries(KratosLayerApplication b64)
if(${LAYER_APP_USE_HDF5} MATCHES TRUE)
    target_link_libraries(KratosLayerApplication ${HDF5_LIBRARIES})
endif()
set_target_properties(KratosLayerApplication PROPERTIES PREFIX "")
install(TARGETS KratosLayerApplication DESTINATION libs )
#######################################################################


## changing the .dll suffix to .pyd for Windows compatibility
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosLayerApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#######################################################################


## install application's python script if required
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
#######################################################################


## Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LayerApplication.py" DESTINATION KratosMultiphysics )
#######################################################################

