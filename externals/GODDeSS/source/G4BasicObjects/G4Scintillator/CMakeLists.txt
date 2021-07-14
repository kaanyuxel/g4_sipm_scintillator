#/*
# * author:      Erik Dietz-Laursonn
# * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
# * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
# */


########## Setup the project ##########
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
# Set project name:
project(G4Scintillator)

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
set(LIBS_GEANT4 ${Geant4_LIBRARIES})

########## Add libraries of the current project ##########
set(LIBS G4PropertyFileReadInTools)

########## The project ##########
# Get list of needed header and source files:
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add include directories to the list of directories searched by the compiler for include files:
include_directories(${PROJECT_SOURCE_DIR}/include)

# Create one library from all source files and add links and dependencies (inside the current project):
add_library(${PROJECT_NAME} ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} ${LIBS_GEANT4} ${LIBS})
add_dependencies(${PROJECT_NAME} ${LIBS})

# add dependencies to the documentation (defined in the top-level cmake file):
add_dependencies(documentation ${PROJECT_NAME})
