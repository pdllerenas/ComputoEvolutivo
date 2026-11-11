find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX)
