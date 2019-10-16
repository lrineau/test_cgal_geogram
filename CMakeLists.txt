project( bench_cgal_geogram )

cmake_minimum_required(VERSION 3.1...3.15)

#Il faut adapter ce path pour qu'il aille chercher l'installation locale.
#Normalement sur windows c'est dans Program Files (x86)\Geogram ou quelque chose de ce goût là
LIST(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/modules")

find_package(Geogram REQUIRED)
find_package(CGAL REQUIRED)

add_executable(_test geogram_bench.cpp )
target_link_libraries(_test PRIVATE CGAL::CGAL Geogram::geogram)

add_executable(test_cgal geogram_bench.cpp )
target_compile_definitions(test_cgal PRIVATE ONLY_CGAL)
target_link_libraries(test_cgal PRIVATE CGAL::CGAL)

add_executable(test_geogram geogram_bench.cpp )
target_compile_definitions(test_geogram PRIVATE ONLY_GEOGRAM)
target_link_libraries(test_geogram PRIVATE Geogram::geogram)

