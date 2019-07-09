#include "spinodal.h"

// Generate array of random numbers between 0 and 1
void rand_ZeroToOne(int nx, int ny, float random_ZeroToOne_array[nx][ny]){
  srand(time(0));
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      random_ZeroToOne_array[i][j] = ( (double)rand() / (double)RAND_MAX );
    }
  }

}
