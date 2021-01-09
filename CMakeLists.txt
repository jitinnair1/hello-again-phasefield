cmake_minimum_required(VERSION 3.14)
project(hello-phasefield)

#Set C Standard
set(CMAKE_C_STANDARD 99)

# Define flags for specific builds
set(CMAKE_C_FLAGS_DEBUG "-Wall -g -O3")

# GCC include directory on macOS which includes openmp.h, openacc.h etc.
include_directories(/usr/local/Cellar/gcc)
include_directories(/usr/local/include)

# Link Libraries
link_directories(/usr/local/lib/)

# Source Files
include_directories(diffusion_eqn_explicit)
include_directories(spinodal_binary_fdm)
include_directories(spinodal_binary_spectral)

add_executable(diffusion_eqn_explicit
        diffusion_eqn_explicit/Diff1DExplicit.c
        diffusion_eqn_explicit/Diff1DExplicit.h
        diffusion_eqn_explicit/write_to_VTK.c)

add_executable(spinodal_binary_fdm
        spinodal_binary_fdm/array_allocate.c
        spinodal_binary_fdm/array_deallocate.c
        spinodal_binary_fdm/free_energy.c
        spinodal_binary_fdm/laplacian.c
        spinodal_binary_fdm/prep_microstructure.c
        spinodal_binary_fdm/rand_ZeroToOne.c
        spinodal_binary_fdm/set_array_zero.c
        spinodal_binary_fdm/solve.c
        spinodal_binary_fdm/solve2.c
        spinodal_binary_fdm/spinodal.c
        spinodal_binary_fdm/spinodal.h
        spinodal_binary_fdm/write_to_VTK.c)

add_executable(spinodal_binary_spectral
        spinodal_binary_spectral/prep_microstructure.c
        spinodal_binary_spectral/rand_ZeroToOne.c
        spinodal_binary_spectral/spinodal_spect.c
        spinodal_binary_spectral/spinodal_spect.h
        spinodal_binary_spectral/write_to_VTK.c
        spinodal_binary_spectral/write_init_conc.c)


# Link Library to Target
target_link_libraries(diffusion_eqn_explicit m)
target_link_libraries(spinodal_binary_fdm m)
target_link_libraries(spinodal_binary_spectral fftw3 m gsl gslcblas)





