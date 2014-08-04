#include "black_scholes.h"
#include "gaussian.h"
#include "random.h" 
#include "util.h"
#include "timer.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>


/**
 * This function is what you compute for each iteration of
 * Black-Scholes.  You don't have to understand it; just call it.
 * "gaussian_random_number" is the current random number (from a
 * Gaussian distribution, which in our case comes from gaussrand1()).
 */
static inline double
black_scholes_value(const double S,
        const double E,
        const double r,
        const double sigma,
        const double T,
        const double gaussian_random_number) {
    const double current_value = S * exp((r - (sigma * sigma) / 2.0) * T +
            sigma * sqrt(T) * gaussian_random_number);
    return exp(-r * T) *
            ((current_value - E < 0.0) ? 0.0 : current_value - E);
    /* return exp (-r * T) * max_double (current_value - E, 0.0); */
}

/**
 * Compute the standard deviation of trials[0 .. M-1].
 */
static double
black_scholes_stddev(int m, double mean, double* trials) {

    const int M = m;
    double variance = 0.0;
    int k;

    for (k = 0; k < M; k++) {
        const double diff = trials[k] - mean;
        /*
         * Just like when computing the mean, we scale each term of this
         * sum in order to avoid overflow.
         */
        variance += diff * diff / (double) M;
    }

    return sqrt(variance);
}

/**
 * Take a pointer to a black_scholes_args_t struct, and return NULL.
 * (The return value is irrelevant, because all the interesting
 * information is written to the input struct.)  This function runs
 * Black-Scholes iterations, and computes the local part of the mean.
 */
static double black_scholes_iterate(int M, double* trials) {

    /* Unpack the IN/OUT struct */

    /* IN (read-only) parameters */
   
    double mean = 0.0;

    /* Temporary variables */
    int k;

    /* Do the Black-Scholes iterations */
    for (k = 0; k < M; k++) {
        /*
         * We scale each term of the sum in order to avoid overflow. 
         * This ensures that mean is never larger than the max
         * element of trials[0 .. M-1].
         */
        mean += trials[k] / (double) M;
    }

    /* 
     * We do the standard deviation computation as a second operation.
     */
    return mean;
}

void *fillTrials(void* trailDetails) {
    
    black_scholes_args_parallel_t* args = (black_scholes_args_parallel_t*) trailDetails;
    gaussrand_state_t gaussrand_state;
    void* prng_stream = NULL;

    const int S = args->S;
    const int E = args->E;
    const double r = args->r;
    const double sigma = args->sigma;
    const double T = args->T;

    const int threadId = args->threadId;
    const int index = args->index;
    int quota = args->quota;
    /* Spawn a random number generator */
    init_prng (index);
    prng_stream = spawn_prng_stream(index);

    /* Initialize the Gaussian random number module for this thread */
    init_gaussrand_state(&gaussrand_state);
    //srand(gettimeofday(NULL, NULL));

    /* Do the Black-Scholes iterations */
    int seed = index *1000;
    int k;
    for (k = 0 * quota; k < quota; k++) {
//        const double gaussian_random_number = gaussrand1(&uniform_random_double,
//                prng_stream,
//                &gaussrand_state);
        const double gaussian_random_number = ((double)rand_r(&seed)/(double)RAND_MAX);
        args->trials[k] = black_scholes_value(S, E, r, sigma, T,
                gaussian_random_number);
    }
}

void
black_scholes_parallel_2(confidence_interval_t* interval,
        const double S,
        const double E,
        const double r,
        const double sigma,
        const double T,
        const int M,
        const int nThreads) {
    pthread_t meanCalThreads[nThreads];
    double mean = 0.0;
    double stddev = 0.0;
    double conf_width = 0.0;
    double variance = 0.0;
    assert(M > 0);
    int quota = (M /nThreads) + 1;
    int i = 0;
    double* trials = (double*) malloc((quota)*(nThreads) * sizeof (double)); 
    black_scholes_args_parallel_t argsArr[nThreads];
    for (i = 0; i < nThreads; i++) {
        double* trialsPerThread = (double*) malloc((quota) * sizeof (double));
        argsArr[i].S = S;
        argsArr[i].E = E;
        argsArr[i].r = r;
        argsArr[i].sigma = sigma;
        argsArr[i].T = T;
        argsArr[i].M = M;
        argsArr[i].mean = 0.0;
        argsArr[i].variance = &variance;
        argsArr[i].index = i;
        argsArr[i].quota = quota;
        argsArr[i].threadId = i;
        argsArr[i].globalMean = &mean;
        argsArr[i].trials = trialsPerThread;
    }

    int k = 0;
    for (k = 0; k < nThreads; k++) {
        int rc = pthread_create(&meanCalThreads[k], NULL, fillTrials, (void *) &argsArr[k]);
        if (rc) {
            printf("Error occurred and the error code is %d\n", rc);
        }
    }

    int j = 0;
    for (j = 0; j < nThreads; j++) {
        pthread_join(meanCalThreads[j], NULL);
    }
    for(j=0; j<nThreads; j++){
       memcpy(trials + (quota)*j, argsArr[j].trials, quota* sizeof (double));   
    }

    mean = black_scholes_iterate(M, trials);
    stddev = black_scholes_stddev(M, mean, trials);

    conf_width = 1.96 * stddev / sqrt((double) M);
    interval->min = mean - conf_width;
    interval->max = mean + conf_width;
}

void deinit_trials(double * trials) {
    if (trials != NULL) {
        free(trials);
    }
}



