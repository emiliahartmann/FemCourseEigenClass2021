# @file FEMCourseEigen/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 3.11.0)

project (FemCourseRef)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#[[
If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
to Debug prior to calling PROJECT(). This assumes you're only build one type of solution
a time. If you're working with multiple solutions like Visual Studio Professional can do,
you need to comment this out.
]]
if(DEFINED CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: Debug, Release.")
else() # CMAKE_BUILD_TYPE
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: Debug, Release.")
endif() # CMAKE_BUILD_TYPE
 
set(CMAKE_CXX_STANDARD 17)

# Eigen
include(cmake/EnableEigen.cmake)
enable_eigen()

# Catch2
include(cmake/EnableCatch2.cmake)
enable_catch2()



## Some definitions for windows
if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	# disabling ALL warnings
	# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
	# disabling SOME warnings
	if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
		set(MSVC_DISABLED_WARNINGS_LIST
			"C4018" # 'expression' : signed/unsigned mismatch
			"C4057" # 'operator' : 'identifier1' indirection to
			# slightly different base types from 'identifier2'
			"C4100" # 'identifier' : unreferenced formal parameter
			"C4127" # conditional expression is constant
			"C4146" # unary minus operator applied to unsigned type,
			# result still unsigned
			"C4244" # 'argument' : conversion from 'type1' to 'type2',
			# possible loss of data
			"C4245" # 'conversion' : conversion from 'type1' to 'type2',
			# signed/unsigned mismatch
			"C4267" # 'var' : conversion from 'size_t' to 'type',
			# possible loss of data
			"C4389" # 'operator' : signed/unsigned mismatch
			"C4706" # assignment within conditional expression
			"C4996" # The POSIX name for this item is deprecated.
			# Instead, use the ISO C and C++ conformant name
		)
		string(REPLACE "C" " -wd" MSVC_DISABLED_WARNINGS_STR ${MSVC_DISABLED_WARNINGS_LIST})
		string(REGEX REPLACE "[/-]W[1234][ ]?" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -MP -W3 ${MSVC_DISABLED_WARNINGS_STR}") # To disable ALL warnings, you can change it to -W0
	endif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
	
endif()


# LINKING FAD
option(USING_FAD "Whether the FAD package will be linked with the library" OFF)

if(USING_FAD)
	add_subdirectory("Fad")
endif(USING_FAD)









## PROJECT FILES
#add_subdirectory(headers)
add_subdirectory(headers)
#add_subdirectory("src vazio")
add_subdirectory(sources)
add_subdirectory(tinyxml)


#Specify where the headers files must be found
include_directories(
	"${PROJECT_SOURCE_DIR}/headers"
	"${PROJECT_SOURCE_DIR}/tinyxml"
)

## inserting main definitions
#specifying the path to FemSC source code
#set(FemSCSOURCE ${FemSC_SRCS} ${FemSC_tinyxml} ${FemSC_HEAD})

set(FemSOURCE ${BASIC_TYPES} ${FemSC_tinyxml} ${HEADERS})

add_subdirectory(mainprograms)
