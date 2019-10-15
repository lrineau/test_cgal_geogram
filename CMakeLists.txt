project( Function_detection_test )

cmake_minimum_required(VERSION 3.1)

#Il faut adapter ce path pour qu'il aille chercher l'installation locale.
#Normalement sur windows c'est dans Program Files (x86)\Geogram ou quelque chose de ce goût là
LIST(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/modules")

find_package(Geogram)
if ( GEOGRAM_FOUND )
  message("GEOGRAM FOUND")
else()
  message("NEEDS GEOGRAM")
endif()

find_package(CGAL REQUIRED)
if(NOT CGAL_FOUND)
  message("NEEDS CGAL")
else()
  add_executable(_test geogram_bench.cpp )
  target_link_libraries(_test CGAL::CGAL Geogram::geogram)
endif()

