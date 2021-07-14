#/*
# * author:      Erik Dietz-Laursonn
# * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
# * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
# */


########## Setup the project ##########
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
# Set project name:
project(G4BasicObjects)

########## The project ##########
# Traverse source folder and find extra includes.
file(GLOB_RECURSE files ABSOLUTE ${CMAKE_CURRENT_SOURCE_DIR} *.h *.hh)

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

# Traverse source folder again and set subprojects.
file(GLOB_RECURSE files ABSOLUTE ${CMAKE_CURRENT_SOURCE_DIR} CMakeLists.txt)

# Add subdirectories which contain CMakeLists, too (they will be executed after this project or, if this project depends on them, before this project):
foreach(file ${files} )
	get_filename_component(path ${file} PATH)

	# Check if the path is the current directory.
	string(COMPARE EQUAL "${path}" "${CMAKE_CURRENT_SOURCE_DIR}" is_current_source_dir)

	if(NOT ${is_current_source_dir})
		# Add cmake subdirectory if it exits.
		if(IS_DIRECTORY ${path})
			message(STATUS "Added subproject: ${path}")
			add_subdirectory(${path})
		endif()
	endif()
endforeach()
