cmake_minimum_required(VERSION 2.8)

project(comet)

# ###################################################################### #
# ###### Set up Rmath dependencies ##################################### #
# ###################################################################### #

set(R_MATH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../extern/Rmath")
set(R_MATH_HEADERS "${R_MATH_DIR}/Rmath")
set(R_MATH_EXT "${R_MATH_DIR}/Rmath/R_ext")

set(R_MATH_LIB "${R_MATH_DIR}/libRmath.a")

message("R_MATH_DIR=${R_MATH_DIR}")
message("R_MATH_HEADERS=${R_MATH_HEADERS}")
message("R_MATH_EXT=${R_MATH_EXT}")

message("R_MATH_LIB=${R_MATH_LIB}")

include_directories(${R_MATH_DIR})
include_directories(${R_MATH_HEADERS})
include_directories(${R_MATH_EXT})

FILE(GLOB RMATH_CFiles "${R_MATH_DIR}/*.c")
list(REMOVE_ITEM RMATH_CFiles "${R_MATH_DIR}/randmtzig.c")
set(CMAKE_C_FLAGS "--std=gnu99 ${CMAKE_C_FLAGS}")

# ###################################################################### #
# ###### Include statements ############################################ #
# ###################################################################### #

include_directories(CometSearch)
include_directories(CometSearch/PeptideFragmentation)
include_directories(CometSearch/PeptideFragmentation/MakeMS2)
include_directories(CometSearch/PeptideFragmentation/MakeMS2/FragmentModels)
include_directories(MSToolkit/include)
include_directories(MSToolkit/src/MSToolkit)
include_directories(MSToolkit/src/expat-2.2.0)
include_directories(MSToolkit/src/mzParser)
include_directories(MSToolkit/src/zlib-1.2.8)

if (APPLE)
        set(CMAKE_OSX_SYSROOT "/")
        include_directories("/usr/local/include/eigen3/")
else ()
	include_directories("/usr/include/eigen3")
endif ()

# ###################################################################### #
# ###### Define Sources ################################################ #
# ###################################################################### #

aux_source_directory(. SOURCE_DIR)
aux_source_directory(CometSearch COMETSEARCH_DIR)
aux_source_directory(CometSearch/PeptideFragmentation PEPTIDEFRAGMENTATION_DIR)
aux_source_directory(CometSearch/PeptideFragmentation/MakeMS2 MAKEMS2_DIR)
aux_source_directory(CometSearch/PeptideFragmentation/MakeMS2/FragmentModels/ MAKEMS2_DIR)
aux_source_directory(MSToolkit/src/expat-2.2.0 EXPAT_DIR)
aux_source_directory(MSToolkit/src/MSToolkit MSTOOKLIT_DIR)
aux_source_directory(MSToolkit/src/mzParser MZPARSER_DIR)
aux_source_directory(MSToolkit/src/zlib-1.2.8 ZLIB_DIR)
aux_source_directory(${R_MATH_DIR} R_MATH_DIR)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -pthread -D__LINUX__ -DCURL_STATICLIB -DHTTP_ONLY" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_NO_THERMORAW -D_NOSQLITE" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include -I/usr/local/include/eigen3 -I/usr/local/include/boost")

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpermissive -pthread -D_NOSQLITE -DHAVE_EXPAT_CONFIG_H -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC")

# Added by Phil 2019-07-15
if (APPLE)
        add_link_options("-v")
        set(CMAKE_CXX_STANDARD 11)

        # 2022-08-02 (Phil):
        # Additional linker flags necessary to get comet to build on my local machine [Mac OS Monterey 12.3.1]
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/opt/llvm/lib/ -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/")
else ()
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -static -lpthread -lm -ldl ")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -lpthread -lm -ldl ")
endif ()

add_executable(comet.exe ${SOURCE_DIR} ${COMETSEARCH_DIR} ${MAKEMS2_DIR} ${PEPTIDEFRAGMENTATION_DIR} ${EXPAT_DIR} ${MSTOOKLIT_DIR} ${MZPARSER_DIR} ${ZLIB_DIR})

add_library(R_MATH_LIB STATIC ${RMATH_CFiles})
set_target_properties(R_MATH_LIB PROPERTIES LINKER_LANGUAGE C)

# ###################################################################### #
# ###### Link libraries ################################################ #
# ###################################################################### #

target_link_libraries(comet.exe ${R_MATH_LIB})
