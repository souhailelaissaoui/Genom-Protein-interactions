cmake_minimum_required(VERSION 2.6)

project(genom-2 CXX)

option(main "Build main version." ON)
option(test "Build tests." ON) #declares a boolean that tells that we are doing a test version.
option(Qt "Building Graphic Interface." ON)

set(CMAKE_BUILD_TYPE Debug) #Allows to debug
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -pedantic -O3 -std=c++11")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if (main)
include_directories(src)
link_directories("/usr/local/lib") #Tells cmake where some libs could potentially be located
link_directories("/usr/lib")
add_executable(Main src/genomic_coordinates.cpp src/main.cpp src/Matrix.cpp src/procedures.cpp src/user_interaction.cpp src/utility.cpp)

target_link_libraries(Main png pthread)

endif (main)

if (test)
	enable_testing()
	find_package(GTest REQUIRED)
	include_directories(${GTEST_INCLUDE_DIRS})
	add_executable(gtests src/genomic_coordinates.cpp test/gTest.cpp src/Matrix.cpp src/user_interaction.cpp src/utility.cpp )
	target_link_libraries(gtests ${GTEST_BOTH_LIBRARIES} pthread) 
	add_test(genom-2 gtests)
endif (test)

find_package(Qt4)
if (QT4_FOUND AND Qt)
	set(CMAKE_AUTOMOC ON)
	set(CMAKE_AUTOOUIC ON)
	execute_process(COMMAND bash ../src/qt.sh )
	include_directories(src /usr/include/qt4/Qt)
	link_directories("/usr/local/lib" "/usr/lib") #Tells cmake where some libs could potentially be located
	add_executable(genom-2 src/askBaseProb.cpp src/genomic_coordinates.cpp src/mainQt.cpp src/Matrix.cpp src/matrixFromSequence.cpp src/procedures.cpp src/resultsWindow.cpp src/sequenceFromMatrix.cpp src/user_interactionQt.cpp src/utility.cpp src/window.cpp)
	target_link_libraries(genom-2 Qt4::QtGui png pthread)
endif (QT4_FOUND AND Qt)

find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)



#To execute CMake: 
#- change "My_first_project" with the project's name
#- change "PopulationSimul" with the program you wan't to compile with CMake 
#- change "src/genome.cpp .... with all the files Cmake will need to compile your program
