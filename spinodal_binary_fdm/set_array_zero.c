#include "spinodal.h"

void set_array_zero(int nx, int ny, float input_array[nx][ny]){
  for (int i = 0; i <= nx; i++) {
    for (int j = 0; j <= ny; j++) {
      input_array[i][j]=0.0;
    }
  }
}
