// dot-product.c
// compile as: gcc dot-product.c [-DNUM_TRIALS=x] -O3 -lrt -lpthread -o prog

#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>		// for clock_gettime()
#include <errno.h>		// for perror()

#define MAX_THREADS     32

#define GET_TIME(x);	if (clock_gettime(CLOCK_MONOTONIC, &(x)) < 0) \
				{ perror("clock_gettime( ):"); exit(EXIT_FAILURE); }
#ifndef NUM_TRIALS
#define NUM_TRIALS  (100*1000*1000)  // problem size
#endif

double  *vec1, *vec2;
double  dot_product_sequential = 0, dot_product_parallel = 0;
long    num_threads = 0;
pthread_mutex_t mutex;

void    initialize_vector_arrays(void);
void    execute_sequential_version(void);
void    handle_sequential_version(void);
void    *run(void *);
void    execute_parallel_version(void);
void    handle_parallel_version(void);
void    handle_parallel_version_with_verification(void);
float   elapsed_time_msec(struct timespec *, struct timespec *, long *, long *);
void    validate_thread_count(char *thread_count);

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
    free(vec1);
    free(vec2);
	return 0;
}

//-----------------sequential code related------------------

void handle_sequential_version()
{
    struct timespec 	start_time_calculation, start_time_total, end_time;
    float 			    comp_time_calculation, comp_time_total; 	// in milli seconds
    unsigned long 		sec, nsec;

    GET_TIME(start_time_total);
	initialize_vector_arrays();

	GET_TIME(start_time_calculation);
	execute_sequential_version();
    GET_TIME(end_time);

    comp_time_calculation = elapsed_time_msec(&start_time_calculation, &end_time, &sec, &nsec);
    comp_time_total = elapsed_time_msec(&start_time_total, &end_time, &sec, &nsec);

    printf("Dot-Product-Sequential: %f : #Trials=%ld : Elapsed-time-for-calculation-only(ms)=%.2f : Elapsed-time-total(ms)=%.2f\n", dot_product_sequential, (long) NUM_TRIALS, comp_time_calculation, comp_time_total);
}

void execute_sequential_version(void)
{
    dot_product_sequential = 0;
    long    i;
    for (i=0; i < NUM_TRIALS; i++) {
        dot_product_sequential += vec1[i] * vec2[i];
    }
}

//-----------------parallel code related------------------

void handle_parallel_version()
{
    struct timespec 	start_time_calculation, start_time_total, end_time;
    float 			    comp_time_calculation, comp_time_total; 	// in milli seconds
    unsigned long 		sec, nsec;

    GET_TIME(start_time_total);
	initialize_vector_arrays();

	GET_TIME(start_time_calculation);
	execute_parallel_version();
    GET_TIME(end_time);

    comp_time_calculation = elapsed_time_msec(&start_time_calculation, &end_time, &sec, &nsec);
    comp_time_total = elapsed_time_msec(&start_time_total, &end_time, &sec, &nsec);

    printf("Dot-Product-Parallel: %f : #Trials=%ld : #Threads=%ld : Elapsed-time-for-calculation-only(ms)=%.2f : Elapsed-time-total(ms)=%.2f\n", dot_product_parallel, (long) NUM_TRIALS, num_threads, comp_time_calculation, comp_time_total);
}

void handle_parallel_version_with_verification()
{
    struct timespec     start_time_vector_init, end_time_vector_init;
    struct timespec 	start_time_calculation_sequntial, end_time_calculation_sequntial;
    struct timespec     start_time_calculation_parallel, end_time_calculation_parallel;
    float 			    comp_time_calculation_sequntial, comp_time_calculation_parallel, comp_time_vector_init; 	// in milli seconds
    unsigned long 		sec, nsec;

    GET_TIME(start_time_vector_init);
	initialize_vector_arrays();
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

    printf("#Trials=%ld : Elapsed-time-for-vector-initialization(ms)=%.2f : Elapsed-time-for-sequential-calculation(ms)=%.2f : Elapsed-time-for-parallel-calculation(ms)=%.2f : #Threads=%ld\n", (long) NUM_TRIALS, comp_time_vector_init, comp_time_calculation_sequntial, comp_time_calculation_parallel, num_threads);
    printf("Dot-Product-Sequential: %f : Dot-Product-Parallel: %f : Difference(Sequential-Parallel): %f\n", dot_product_sequential, dot_product_parallel, (dot_product_sequential - dot_product_parallel));
}

void execute_parallel_version()
{
    pthread_t 		    tid[MAX_THREADS];
    int 			    i, myid[MAX_THREADS];
    dot_product_parallel = 0;

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
    long	i;
    double  my_dot_product=0;

	for (i=start; i < end; i++) {
		my_dot_product += vec1[i] * vec2[i];
	}

   	pthread_mutex_lock (&mutex);
   	dot_product_parallel += my_dot_product;
   	pthread_mutex_unlock (&mutex);

   	pthread_exit(NULL);
}


//-----------------utility functions------------------

void validate_thread_count(char *thread_count)
{
    num_threads = atoi(thread_count);
    if (num_threads <= 0 || num_threads > MAX_THREADS) {
        printf("num_threads should be in range 1 and %d(MAXTHREADS)\n", MAX_THREADS);
        exit(0);
    }
}

void initialize_vector_arrays(void)
{
	long i;

    vec1 = (double *)malloc(sizeof(double)*NUM_TRIALS);
    vec2 = (double *)malloc(sizeof(double)*NUM_TRIALS);

    for (i=0; i < NUM_TRIALS; i++) {
        // adding 1 because we need number between 1 and 2.
        vec1[i] = (((double) rand())/ RAND_MAX) + 1;
    	vec2[i] = (((double) rand())/ RAND_MAX) + 1;
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
