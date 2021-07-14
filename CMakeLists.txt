cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(sipm_scintillator_app)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# Boost
set(Boost_USE_MULTITHREADED OFF)
find_package(Boost 1.50.0 COMPONENTS date_time program_options filesystem system regex REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
set(LIBS ${LIBS} ${Boost_LIBRARIES})

# Zlib 
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ZLIB_LIBRARIES})

add_subdirectory(externals/GODDeSS/source)
include_directories(externals/GODDeSS/source/G4BasicObjects/G4PropertyFileReadInTools/include)
include_directories(externals/GODDeSS/source/G4BasicObjects/G4Scintillator/include)
include_directories(externals/GODDeSS/source/G4BasicObjects/G4Fibre/include)
include_directories(externals/GODDeSS/source/G4BasicObjects/G4GeneralParticleSource/include)
include_directories(externals/GODDeSS/source/G4BasicObjects/G4PhotonDetector/include)
include_directories(externals/GODDeSS/source/G4BasicObjects/DataProcessing/include)



add_subdirectory(externals/g4sipm)
include_directories(externals/g4sipm/g4sipm/include)
include_directories(externals/g4sipm/externals/jansson/src)
include_directories(externals/g4sipm/externals/gtest/include)

# We presume the existence of four directories
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB_RECURSE macros RELATIVE ${PROJECT_SOURCE_DIR} macros/*.mac)
file(GLOB_RECURSE materialproperties RELATIVE ${PROJECT_SOURCE_DIR} externals/GODDeSS/MaterialProperties/*)


# Enable macros for out-of-source build
foreach(_file ${macros})
  configure_file(
    ${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
    )
endforeach()
# Enable Material Properties for out-of-source build
foreach(_file ${materialproperties})
  configure_file(
    ${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
    )
endforeach()
#Copy scripts
configure_file( 
    ${PROJECT_SOURCE_DIR}/DrawHist.cpp
    ${PROJECT_BINARY_DIR}/DrawHist.cpp 
    COPYONLY
    )


set(LIBS ${LIBS} boost_program_options g4sipm boost_date_time janssonset G4PropertyFileReadInTools G4Scintillator G4Fibre G4GeneralParticleSource G4PhotonDetector DataProcessing)
add_executable(sipm_scintillator_app main.cc ${sources} ${headers})
target_link_libraries(sipm_scintillator_app ${Geant4_LIBRARIES} g4sipm boost_date_time jansson G4PropertyFileReadInTools G4Scintillator G4Fibre G4GeneralParticleSource G4PhotonDetector DataProcessing)

install(TARGETS sipm_scintillator_app DESTINATION bin)
