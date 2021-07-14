#/*
# * author:      Erik Dietz-Laursonn
# * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
# * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
# */


########## Setup the project ##########
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
# Set project name:
project(RunSimulation)

########## Add external libraries ##########
#------ Geant4 ------#
# Find Geant4 package, activating all available UI and Vis drivers by default
option(WITH_GEANT4_UIVIS "Build with all Geant4 UIs and all visualization drivers" ON)
option(WITH_GEANT4_OPENGL "Build with Geant4 shell UI and OpenGl visualization" OFF)
option(WITH_GEANT4_QT "Build with Geant4 Qt UI" OFF)

if(WITH_GEANT4_OPENGL)
	find_package(Geant4 REQUIRED vis_opengl_x11 ui_tcsh)
elseif(WITH_GEANT4_QT)
	find_package(Geant4 REQUIRED qt)
elseif(WITH_GEANT4_UIVIS)
	find_package(Geant4 REQUIRED ui_all vis_all)
else()
	find_package(Geant4 REQUIRED)
endif()

# Add Geant4 include directories and compile definitions to the list of directories searched by the compiler for include files:
include(${Geant4_USE_FILE})

# Set Geant4 library list:
# set(LIBS_GEANT4 ${Geant4_LIBRARIES})   // NOTE: These libraries are already included by the G4BasicObjects of GODDeSS and an additional inclusion here leads to linking errors since Geant4.10!

#------ Boost ------#
# Find Boost package
find_package(Boost REQUIRED regex)

# Add include directories to the list of directories searched by the compiler for include files:
include_directories(${Boost_INCLUDE_DIRS})

# Set Boost library list
set(LIBS_BOOST boost_regex)

#------ zlib ------#
# Find zlib package
find_package(ZLIB REQUIRED)

# Add include directories to the list of directories searched by the compiler for include files:
include_directories(${ZLIB_INCLUDE_DIRS})

# Set zlib library list
set(LIBS_ZLIB ${ZLIB_LIBRARIES})

########## Add libraries of the current project ##########
set(LIBS G4PropertyFileReadInTools G4Scintillator G4Fibre G4GeneralParticleSource G4PhotonDetector DataProcessing)

########## The project ##########
# Get list of needed header and source files:
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/DataProcessing/*.cc ${PROJECT_SOURCE_DIR}/src/Preparation/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/DataProcessing/*.hh ${PROJECT_SOURCE_DIR}/include/Preparation/*.hh)

# Add include directories to the list of directories searched by the compiler for include files:
include_directories(${PROJECT_SOURCE_DIR}/include/DataProcessing ${PROJECT_SOURCE_DIR}/include/Preparation)

# Create the executable and add links and dependencies (inside the current project):
add_executable(${PROJECT_NAME} RunSimulation.cc ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${LIBS_GEANT4} ${LIBS_BOOST} ${LIBS_ZLIB} ${LIBS})
add_dependencies(${PROJECT_NAME} ${LIBS})

# Install the executable to the destination directory:
set(INSTALL_DIRECTORY ${CMAKE_BINARY_DIR})
install(TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_DIRECTORY})

# add dependencies to the documentation (defined in the top-level cmake file):
add_dependencies(documentation ${PROJECT_NAME})
