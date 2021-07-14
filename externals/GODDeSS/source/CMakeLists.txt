#/*
# * author:      Erik Dietz-Laursonn
# * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
# * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
# */


########## Setup the project ##########
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
# Set project name:
project(Simulation)

########## The project ##########
# Traverse source folder and find extra includes.
file(GLOB_RECURSE files ABSOLUTE ${CMAKE_CURRENT_SOURCE_DIR} *.h *.hh)

option(BUILD_SHARED_LIBS "Build package with shared libraries." OFF)

set(last_include_path "")

foreach(file ${files} )
	get_filename_component(path ${file} PATH)

	# Check if the path has not been processed before.
	string(COMPARE EQUAL "${path}" "${last_include_path}" already_included)
	if(NOT ${already_included})
		# Adapt last_include_path.
		set(last_include_path "${path}")
		# Add include path if include subdirectory exits.
		if(IS_DIRECTORY ${path})
			message(STATUS "Added include: ${path}")
			include_directories(${path})
		endif()
	endif()
endforeach()

# make documentation (after all subproject (which add dependencies to the documentation) are build):
add_custom_target(documentation ALL
				COMMAND echo Creating documentation
				COMMAND rm -rf ../documentation
				COMMAND mkdir -p ../documentation
				COMMAND doxygen GODDeSS.dox > /dev/null
				WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
# add_custom_command(TARGET documentation
#                      POST_BUILD
#                      COMMAND echo ARGS Creating documentation pdf-files
#                      COMMAND sed ARGS -e "s#makeindex#makeindex -q#" -i Makefile
#                      COMMAND make all > /dev/null
#                      COMMAND cp refman.pdf ../documentation.pdf
#                      COMMAND touch ../documentation.html
#                      COMMAND echo ARGS \"<html> <head> <meta http-equiv=\\\"refresh\\\" content=\\\"0\\; URL=html/index.html\\\"/> </head> </html>\" > ../documentation.html
#                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../documentation/latex)

# Add subdirectories which contain CMakeLists, too (they will be executed after this project or, if this project depends on them, before this project):
message(STATUS "Added subproject: ${CMAKE_CURRENT_SOURCE_DIR}/G4BasicObjects")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/G4BasicObjects)
