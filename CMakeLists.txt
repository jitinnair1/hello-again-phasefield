cmake_minimum_required(VERSION 3.14)
project(phasefieldC)

#Set C Standard
set(CMAKE_C_STANDARD 99)

# Define flags for specific builds
set(CMAKE_C_FLAGS_DEBUG "-Wall -g -O2")

# GCC include directory on macOS which includes openmp.h, openacc.h etc.
include_directories(/opt/local/lib/gcc9/gcc/x86_64-apple-darwin18/9.1.0/include/)
include_directories(/usr/local/include)

# Link Libraries
link_directories(/usr/local/lib/)

# Source Files
include_directories(1d_diff_explicit)
include_directories(spinodal)
include_directories(spinodal_spect)
include_directories(spinodal_ternary)

add_executable(1d_diff_explicit
        1d_diff_explicit/Diff1DExplicit.c
        1d_diff_explicit/Diff1DExplicit.h
        1d_diff_explicit/write_to_VTK.c)

add_executable(spinodal
        spinodal/array_allocate.c
        spinodal/array_deallocate.c
        spinodal/free_energy.c
        spinodal/laplacian.c
        spinodal/prep_microstructure.c
        spinodal/rand_ZeroToOne.c
        spinodal/set_array_zero.c
        spinodal/solve.c
        spinodal/solve2.c
        spinodal/spinodal.c
        spinodal/spinodal.h
        spinodal/write_to_VTK.c)

add_executable(spinodal_spect
        spinodal_spect/array_allocate.c
        spinodal_spect/array_deallocate.c
        spinodal_spect/prep_microstructure.c
        spinodal_spect/rand_ZeroToOne.c
        spinodal_spect/spinodal_spect.c
        spinodal_spect/spinodal_spect.h
        spinodal_spect/write_to_VTK.c)

add_executable(spinodal_ternary
        spinodal_ternary/array_allocate.c
        spinodal_ternary/array_deallocate.c
        spinodal_ternary/prep_microstructure.c
        spinodal_ternary/rand_ZeroToOne.c
        spinodal_ternary/ternary_coarsening.c
        spinodal_ternary/ternary_coarsening.h
        spinodal_ternary/write_to_VTK.c)

# Link Library to Target
target_link_libraries(1d_diff_explicit m)
target_link_libraries(spinodal m)
target_link_libraries(spinodal_spect fftw3 m)
target_link_libraries(spinodal_ternary fftw3 m)





