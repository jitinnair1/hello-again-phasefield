#include "spinodal.h"

void array_deallocate(int ny, float **array) {
  for (int i = 0; i < ny; i++)
		free(array[i]);
	free(array);
}
