################################################################################
#   CMakeLists.txt for wmc4                                                    #
#                                                                              #
#                      - 20. Oct. 2018. Hoyong Jeong (hyjeong@hep.korea.ac.kr) #
################################################################################

#------------------------------------------------------------------------------#
#   Setup the project                                                          #
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(WMC4)

#------------------------------------------------------------------------------#
#   Find Geant4 package, activating all available UI and Vis drivers by        #
# default. You can set WITH_GEANT4_UIVIS to OFF via the command line or        #
# ccmake/cmake-gui to build a batch mode only executable.                      #
#------------------------------------------------------------------------------#
option(WITH_GEANT4_UIVIS "Build wmc4 with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
	find_package(Geant4 REQUIRED ui_all vis_all)
else()
	find_package(Geant4 REQUIRED)
endif()

#------------------------------------------------------------------------------#
#   Find ROOT package.                                                         #
#------------------------------------------------------------------------------#
find_package(ROOT REQUIRED EG)

#------------------------------------------------------------------------------#
#   Setup Geant4 include directories and compile definitions                   #
#   Setup include directory for this project                                   #
#------------------------------------------------------------------------------#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#------------------------------------------------------------------------------#
#   Setup ROOT include directories and compile definitions.                    #
#------------------------------------------------------------------------------#
include(${ROOT_USE_FILE})

#------------------------------------------------------------------------------#
#   Setup PLUTO using user's environment variable "PLUTOSYS"                   #
# You can edit directory path depending on your system.                        #
#------------------------------------------------------------------------------#
link_directories($ENV{PLUTOSYS})
include_directories($ENV{PLUTOSYS}/include)

#------------------------------------------------------------------------------#
#   Locate sources and headers for this project                                #
#   NB: headers are included so they will show up in IDEs                      #
#------------------------------------------------------------------------------#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#------------------------------------------------------------------------------#
#   Add the executable, and link it to the Geant4 libraries                    #
#------------------------------------------------------------------------------#
add_executable(wmc4 wmc4.cc ${sources} ${headers})
target_link_libraries(wmc4 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Pluto)

#------------------------------------------------------------------------------#
#   Copy all scripts to the build directory, i.e. the directory in which we    #
# build wmc4. This is so that we can run the executable directly because it    #
# relies on these scripts being in the current working directory.              #
#------------------------------------------------------------------------------#
set(WMC4_SCRIPTS
	README
	config.cfg
	init_vis.mac
	vis.mac
	wasa_geometry.pdf
)

foreach(_script ${WMC4_SCRIPTS})
	configure_file(
		${PROJECT_SOURCE_DIR}/${_script}
		${PROJECT_BINARY_DIR}/${_script}
		COPYONLY
	)
endforeach()

file(COPY ${PROJECT_SOURCE_DIR}/scripts DESTINATION ${PROJECT_BINARY_DIR})

#------------------------------------------------------------------------------#
#   For internal Geant4 use - but has no effect if you build this wmc4         #
# standalone                                                                   #
#------------------------------------------------------------------------------#
add_custom_target(WMC4 DEPENDS wmc4)

#------------------------------------------------------------------------------#
#   Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX       #
#------------------------------------------------------------------------------#
install(TARGETS wmc4 DESTINATION bin)

#------------------------------------------------------------------------------#
#   Flag setting                                                               #
#------------------------------------------------------------------------------#
set (CMAKE_EXE_LINKER_FLAGS "-Wl,--no-as-needed")
