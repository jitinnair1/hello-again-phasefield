#include "spinodal_spect.h"

fftw_complex** array_fft_allocate(int nx, int ny, float **array){
  array=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
  return array;
}
