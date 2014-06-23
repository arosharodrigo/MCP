// matrix-multiplication.c
// compile as: gcc matrix-multiplication.c [-DNUM_TRIALS=x] -O3 -lrt -lpthread -o prog

#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>		// for clock_gettime()
#include <errno.h>		// for perror()

#define MAX_THREADS     32

#define GET_TIME(x);	if (clock_gettime(CLOCK_MONOTONIC, &(x)) < 0) \
				{ perror("clock_gettime( ):"); exit(EXIT_FAILURE); }
#ifndef NUM_TRIALS
#define NUM_TRIALS  600  // problem size
#endif

double  **mat1, **mat2;
double  **mat_mul_sequential, **mat_mul_parallel;
long    num_threads = 0;
pthread_mutex_t mutex;

void    execute_sequential_version(void);
void    handle_sequential_version(void);
void    *run(void *);
void    execute_parallel_version(void);
void    handle_parallel_version(void);
void    handle_parallel_version_with_verification(void);
float   elapsed_time_msec(struct timespec *, struct timespec *, long *, long *);
void    initialize_matrices(void);
double  calculate_point_value(long x, long y);
double  calculate_cumulative_diff(double **matrix2, double **matrix1);
void    validate_thread_count(char *thread_count);
void    print_matrix(double **matrix);

int main (int argc, char *argv[])
{
    bool sequential_condition = (argc == 2 && (strcmp(argv[1],"-s") == 0));
    bool parallel_condition = (argc == 3 && (strcmp(argv[1],"-p"))  == 0);
    bool verify_condition = (argc == 4 && strcmp(argv[1],"-p")  == 0 && strcmp(argv[3],"-v") == 0);
    srand(time(NULL));

    if (sequential_condition) {
        printf("Sequential version executing. Please wait.....\n");
        handle_sequential_version();
    } else if (parallel_condition) {
        validate_thread_count(argv[2]);
        printf("Parallel version executing. Please wait.....\n");
        handle_parallel_version();
    } else if (verify_condition) {
        validate_thread_count(argv[2]);
        printf("Parallel verification version executing. Please wait.....\n");
        handle_parallel_version_with_verification();
    } else {
        printf("Usage: [%s -s] OR [%s -p <num_threads>] OR [%s -p <num_threads> -v]\n", argv[0], argv[0], argv[0]);
		exit(0);
    }
	return 0;
}

//-----------------sequential code related------------------

void handle_sequential_version()
{
    struct timespec 	start_time_calculation, start_time_total, end_time;
    float 			    comp_time_calculation, comp_time_total; 	// in milli seconds
    unsigned long 		sec, nsec;

    GET_TIME(start_time_total);
	initialize_matrices();

	GET_TIME(start_time_calculation);
	execute_sequential_version();
    GET_TIME(end_time);

    comp_time_calculation = elapsed_time_msec(&start_time_calculation, &end_time, &sec, &nsec);
    comp_time_total = elapsed_time_msec(&start_time_total, &end_time, &sec, &nsec);
    printf("Matrix-Multiplication-Sequential - #Trials=%ld : Elapsed-time-for-calculation-only(ms)=%.2f : Elapsed-time-total(ms)=%.2f\n", (long) NUM_TRIALS, comp_time_calculation, comp_time_total);
}

void execute_sequential_version(void)
{
    long i, j;
    for (i=0; i < NUM_TRIALS; i++) {
        for (j=0; j < NUM_TRIALS; j++) {
            mat_mul_sequential[i][j] = calculate_point_value(i, j);
    	}
    }
}

//-----------------parallel code related------------------

void handle_parallel_version()
{
    struct timespec 	start_time_calculation, start_time_total, end_time;
    float 			    comp_time_calculation, comp_time_total; 	// in milli seconds
    unsigned long 		sec, nsec;

    GET_TIME(start_time_total);
	initialize_matrices();

	GET_TIME(start_time_calculation);
	execute_parallel_version();
    GET_TIME(end_time);

    comp_time_calculation = elapsed_time_msec(&start_time_calculation, &end_time, &sec, &nsec);
    comp_time_total = elapsed_time_msec(&start_time_total, &end_time, &sec, &nsec);

    printf("Matrix-Multiplication-Parallel - #Trials=%ld : #Threads=%ld : Elapsed-time-for-calculation-only(ms)=%.2f : Elapsed-time-total(ms)=%.2f\n", (long) NUM_TRIALS, num_threads, comp_time_calculation, comp_time_total);
}

void handle_parallel_version_with_verification()
{
    struct timespec     start_time_vector_init, end_time_vector_init;
    struct timespec 	start_time_calculation_sequntial, end_time_calculation_sequntial;
    struct timespec     start_time_calculation_parallel, end_time_calculation_parallel;
    float 			    comp_time_calculation_sequntial, comp_time_calculation_parallel, comp_time_vector_init; 	// in milli seconds
    unsigned long 		sec, nsec;
    double              cumulative_diff;

    GET_TIME(start_time_vector_init);
	initialize_matrices();
    GET_TIME(end_time_vector_init);

    GET_TIME(start_time_calculation_sequntial);
    execute_sequential_version();
    GET_TIME(end_time_calculation_sequntial);

	GET_TIME(start_time_calculation_parallel);
	execute_parallel_version();
    GET_TIME(end_time_calculation_parallel);

    comp_time_vector_init = elapsed_time_msec(&start_time_vector_init, &end_time_vector_init, &sec, &nsec);
    comp_time_calculation_sequntial = elapsed_time_msec(&start_time_calculation_sequntial, &end_time_calculation_sequntial, &sec, &nsec);
    comp_time_calculation_parallel = elapsed_time_msec(&start_time_calculation_parallel, &end_time_calculation_parallel, &sec, &nsec);

    cumulative_diff = calculate_cumulative_diff(mat_mul_sequential, mat_mul_parallel);

    printf("#Trials=%ld : Elapsed-time-for-vector-initialization(ms)=%.2f : Elapsed-time-for-sequential-calculation(ms)=%.2f : Elapsed-time-for-parallel-calculation(ms)=%.2f : #Threads=%ld\n", (long) NUM_TRIALS, comp_time_vector_init, comp_time_calculation_sequntial, comp_time_calculation_parallel, num_threads);
    printf("Cumulative-Difference(Parallel-Sequential): %f\n", cumulative_diff);
}

void execute_parallel_version()
{
    pthread_t 		    tid[MAX_THREADS];
    int 			    i, myid[MAX_THREADS];

    pthread_mutex_init(&mutex, NULL);
    for (i = 0; i < num_threads; i++) {
	    myid[i] = i;
		pthread_create(&tid[i], NULL, run, &myid[i]);
    }
    for (i = 0; i < num_threads; i++)  {
	    pthread_join(tid[i], NULL);
    }
	pthread_mutex_destroy(&mutex);
}

void *run(void * tid)
{
    int 	myid = *((int *)tid);
    long 	start = (myid * (long)NUM_TRIALS)/num_threads;
    long 	end = ((myid+1) * (long)NUM_TRIALS) /num_threads;
    long	i, j, size;
    double  **mat_temp;

    size = end - start;
    mat_temp = malloc(size * sizeof(double*));

    for(i = 0; i < size; i++) {
        mat_temp[i] = malloc( NUM_TRIALS * sizeof(double));
    }

	for (i=0; i < size; i++) {
        for (j=0; j < NUM_TRIALS; j++) {
            mat_temp[i][j] = calculate_point_value(i + start, j);
        }
	}

   	pthread_mutex_lock (&mutex);
    for (i=0; i < size; i++) {
        for (j=0; j < NUM_TRIALS; j++) {
            mat_mul_parallel[i + start][j] = mat_temp[i][j];
        }
	}
   	pthread_mutex_unlock (&mutex);
   	pthread_exit(NULL);
}


//-----------------utility functions------------------

// used in both sequential and parallel versions to calculate matrix multiplication of a particular point.
double calculate_point_value(long x, long y)
{
    long k;
    double res = 0;
    for (k=0; k < NUM_TRIALS; k++) {
        res += mat1[x][k] * mat2[k][y];
    }
    return res;
}

double calculate_cumulative_diff(double **matrix1, double **matrix2)
{
    long i, j;
    double res = 0;

    for (i=0; i < NUM_TRIALS; i++) {
        for (j=0; j < NUM_TRIALS; j++) {
    	    res += fabs(matrix2[i][j] - matrix1[i][j]);
    	}
    }
    return res;
}

void validate_thread_count(char *thread_count)
{
    num_threads = atoi(thread_count);
    if (num_threads <= 0 || num_threads > MAX_THREADS) {
        printf("num_threads should be in range 1 and %d(MAXTHREADS)\n", MAX_THREADS);
        exit(0);
    }
}

void print_matrix(double **matrix)
{
    long i, j;
    for (i=0; i < NUM_TRIALS; i++) {
        for (j=0; j < NUM_TRIALS; j++) {
            printf("%f ", matrix[i][j]);
    	}
    	printf("\n");
    }
    printf("\n");
}

void initialize_matrices(void)
{
	long i, j;

    // allocating sizes for all matrices
    mat1 = malloc( NUM_TRIALS * sizeof(double*));
    mat2 = malloc( NUM_TRIALS * sizeof(double*));
    mat_mul_sequential = malloc( NUM_TRIALS * sizeof(double*));
    mat_mul_parallel = malloc( NUM_TRIALS * sizeof(double*));

    for(i = 0; i < NUM_TRIALS; i++) {
        mat1[i] = malloc( NUM_TRIALS * sizeof(double));
        mat2[i] = malloc( NUM_TRIALS * sizeof(double));
        mat_mul_sequential[i] = malloc( NUM_TRIALS * sizeof(double));
        mat_mul_parallel[i] = malloc( NUM_TRIALS * sizeof(double));
    }

    // initialize mat1 and mat2 matrices with random numbers while
    // mat_mul_sequential and mat_mul_parallel initialize to 0
    for (i=0; i < NUM_TRIALS; i++) {
        for (j=0; j < NUM_TRIALS; j++) {
            // adding 1 because we need number between 1 and 2.
            mat1[i][j] = (((double) rand())/ RAND_MAX) + 1;
    	    mat2[i][j] = (((double) rand())/ RAND_MAX) + 1;
    	    mat_mul_sequential[i][j] = 0;
    	    mat_mul_parallel[i][j] = 0;
    	}
    }
}

float elapsed_time_msec(struct timespec *begin, struct timespec *end, long *sec, long *nsec)
{
	if (end->tv_nsec < begin->tv_nsec) {
        *nsec = 1000000000 - (begin->tv_nsec - end->tv_nsec);
    	*sec  = end->tv_sec - begin->tv_sec -1;
    } else {
		*nsec = end->tv_nsec - begin->tv_nsec;
    	*sec  = end->tv_sec - begin->tv_sec;
	}
	return (float) (*sec) * 1000 + ((float) (*nsec)) / 1000000;
}
