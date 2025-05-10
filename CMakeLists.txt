cmake_minimum_required(VERSION 3.25)
project(FR1D)

set(CMAKE_CXX_STANDARD 23)

add_executable(FR1D main.cpp indexing.cpp indexing.h basis.cpp basis.h SailorMath.cpp SailorMath.h EulerFlux.cpp EulerFlux.h SpatialDiscretization.cpp SpatialDiscretization.h initial.cpp initial.h)
