# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(ValidationModule)
find_package(Falaise REQUIRED)

# Build a dynamic library from our sources
add_library(ValidationModule SHARED ValidationModule.h ValidationModule.cpp TrackDetails.h trackDetails.cpp)

# Link it to the FalaiseModule library
# This ensures the correct compiler flags, include paths
# and linker flags are applied to our dynamic library.
# - Explicit link to ROOT Physics library which is not linked
#   by flreconstruct
target_link_libraries(ValidationModule
  PUBLIC
    Falaise::FalaiseModule
    ${ROOT_Physics_LIBRARY}
    )

# Configure example pipeline script for use from the build dir
configure_file("ValidationModuleExample.conf.in" "ValidationModuleExample.conf" @ONLY)

# Add a basic test of reading a brio file output by the
# standard pipeline
enable_testing()
# - Simulate
add_test(NAME testValidationModule_simulate
  COMMAND Falaise::flsimulate -o test-simulate.brio
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
# - Reconstruct
add_test(NAME testValidationModule_reconstruct
  COMMAND Falaise::flreconstruct -i test-simulate.brio -p urn:snemo:demonstrator:reconstruction:1.0.0 -o test-reconstruct.brio
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
set_tests_properties(testValidationModule_reconstruct
  PROPERTIES DEPENDS testValidationModule_simulate
  )
# - Run Module
add_test(NAME testValidationModule_Validation
  COMMAND Falaise::flreconstruct -i test-reconstruct.brio -p ValidationModuleExample.conf
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
set_tests_properties(testValidationModule_Validation
  PROPERTIES DEPENDS testValidationModule_reconstruct
  )
