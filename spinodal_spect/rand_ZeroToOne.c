#include "spinodal_spect.h"

// Generate array of random numbers between 0 and 1
float** rand_ZeroToOne(int nx, int ny, float **random_ZeroToOne_array){
  //srand(time(0));
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      random_ZeroToOne_array[i][j] = ( (double)rand() / (double)RAND_MAX );
    }
  }
  return random_ZeroToOne_array;

}
