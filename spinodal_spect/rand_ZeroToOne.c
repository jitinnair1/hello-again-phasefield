#include "spinodal_spect.h"
#include <gsl/gsl_rng.h>

// Generate array of random numbers between 0 and 1
double *rand_ZeroToOne(int Nx, int Ny, int seed, double *random_ZeroToOne_array) {
        const gsl_rng_type *T;
        gsl_rng *r;

        int i,
        NxNy = Nx * Ny;

        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, seed);

        // Generate array of random numbers between 0 and 1
        for (i = 0; i < NxNy; i++) {
            random_ZeroToOne_array[i] = gsl_rng_uniform(r);
        }

        gsl_rng_free(r);

        return random_ZeroToOne_array;

    }