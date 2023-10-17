list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/EllipticProblem.hpp)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/Fracture3D.hpp)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/P1MatrixAssembler.hpp)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/XfemDofHandler.hpp)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utilities.hpp)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/PhysicalParameters.hpp)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/discontinousTestProblem_1.h)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/testAllZero.hpp)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/testOneZero.hpp)


list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/EllipticProblem.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/Fracture3D.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/P1MatrixAssembler.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/XfemDofHandler.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/Utilities.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/PhysicalParameters.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/discontinousTestProblem_1.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/testAllZero.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/testOneZero.cpp)



list(APPEND xfem_3d_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(xfem_3d_sources ${xfem_3d_sources} PARENT_SCOPE)
set(xfem_3d_headers ${xfem_3d_headers} PARENT_SCOPE)
set(xfem_3d_includes ${xfem_3d_includes} PARENT_SCOPE)
