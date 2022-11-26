cmake_minimum_required(VERSION 3.9)

project(XFEM_3D LANGUAGES CXX C)

# Create variables
###############################################################################################
set(XFEM_3D_SOURCES ${PROJECT_NAME}_sources)
set(XFEM_3D_HEADERS ${PROJECT_NAME}_headers)
set(XFEM_3D_INCLUDE ${PROJECT_NAME}_include)
set(XFEM_3D_LINKED_LIBRARIES ${PROJECT_NAME}_linked_libraries)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

# Add dependencies
###############################################################################################

## GeDiM
find_package(GeDiM CONFIG REQUIRED)
list(APPEND ${XFEM_3D_LINKED_LIBRARIES} PRIVATE GeDiM::GeDiM)

## Blas
find_package(BLAS REQUIRED)

## Lapack
find_package(LAPACK REQUIRED)

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)

## VTK
if(${ENABLE_VTK})
  find_package(VTK CONFIG REQUIRED)
endif()

## triangle
find_package(triangle CONFIG REQUIRED)

## tetgen
find_package(tetgen CONFIG REQUIRED)

# Insert Sources
###############################################################################################
add_subdirectory(src)
list(APPEND ${XFEM_3D_SOURCES} ${xfem_3d_sources})
list(APPEND ${XFEM_3D_HEADERS} ${xfem_3d_headers})
list(APPEND ${XFEM_3D_INCLUDE} ${xfem_3d_includes})

# Create library test executable
###############################################################################################
add_executable(${PROJECT_NAME}
	main.cpp
	${${XFEM_3D_SOURCES}}
    ${${XFEM_3D_HEADERS}})

target_include_directories(${PROJECT_NAME}
    PRIVATE ${${XFEM_3D_INCLUDE}})

target_link_libraries(${PROJECT_NAME} ${${XFEM_3D_LINKED_LIBRARIES}})
