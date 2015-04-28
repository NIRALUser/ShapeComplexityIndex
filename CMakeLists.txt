cmake_minimum_required (VERSION 2.8)
project (Cal_Complexity)

# GenerateCLP
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
generateCLP(CLP_XML Complex.xml)

# MESH library
add_library(Mesh SHARED IMPORTED)

IF(APPLE)
	message("IS THIS Mac OX??")
	set_target_properties(Mesh PROPERTIES IMPORTED_LOCATION ../src/lib/libMesh.1.0.dylib)
ELSE()
	message("IS THIS LINUX or WIN32??")
	set_target_properties(Mesh PROPERTIES IMPORTED_LOCATION ../src/lib/libMesh.so.1.0)
ENDIF(APPLE)	
	
include_directories(../src/include)

# sources
#set(SRC_FILES src/main.cpp)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

CONFIGURE_FILE(Cal_Complex_opt.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Cal_Complex_opt.x.cpp)
CONFIGURE_FILE(surf_reg_model_left.obj template.obj COPYONLY)

# executable file
#add_executable(MeshExample ${SRC_FILES})
add_executable(Cal_Complex Cal_Complex_opt.x.cpp Calculators.h Calculators.cpp transportSimplex.h ${CLP_XML})
target_link_libraries(Cal_Complex Mesh)