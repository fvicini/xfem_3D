list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/EllipticProblem.hpp)
list(APPEND xfem_3d_headers ${CMAKE_CURRENT_SOURCE_DIR}/Pippo.hpp)

list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/EllipticProblem.cpp)
list(APPEND xfem_3d_sources ${CMAKE_CURRENT_SOURCE_DIR}/Pippo.cpp)

list(APPEND xfem_3d_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(xfem_3d_sources ${xfem_3d_sources} PARENT_SCOPE)
set(xfem_3d_headers ${xfem_3d_headers} PARENT_SCOPE)
set(xfem_3d_includes ${xfem_3d_includes} PARENT_SCOPE)
