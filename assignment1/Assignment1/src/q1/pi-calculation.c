// pi-calculation.c
// compile as: gcc pi-calculation.c [-DNUM_TRIALS=x] -O3 -lrt -lpthread -o prog

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

float 	pi_seqential = 0, pi_parallel = 0;
long    num_threads = 0, parallel_hits = 0;
pthread_mutex_t mutex;

void    execute_sequential_version(void);
void    handle_sequential_version(void);
void    *run(void *);
void    execute_parallel_version(void);
void    handle_parallel_version(void);
float   elapsed_time_msec(struct timespec *, struct timespec *, long *, long *);
void    validate_thread_count(char *thread_count);

int main (int argc, char *argv[])
{
    bool sequential_condition = (argc == 2 && (strcmp(argv[1],"-s") == 0));
    bool parallel_condition = (argc == 3 && (strcmp(argv[1],"-p"))  == 0);

    if (sequential_condition) {
        printf("Sequential version executing. Please wait.....\n");
        handle_sequential_version();
    } else if (parallel_condition) {
        validate_thread_count(argv[2]);
        printf("Parallel version executing. Please wait.....\n");
        handle_parallel_version();
    } else {
        printf("Usage: [%s -s] OR [%s -p <num_threads>]\n", argv[0], argv[0], argv[0]);
		exit(0);
    }
	return 0;
}

//-----------------sequential code related------------------

void handle_sequential_version()
{
    struct timespec 	start_time, end_time;
    float 			    comp_time_total; 	// in milli seconds
    unsigned long 		sec, nsec;

    GET_TIME(start_time);
	execute_sequential_version();
    GET_TIME(end_time);

    comp_time_total = elapsed_time_msec(&start_time, &end_time, &sec, &nsec);

    printf("Pi-Sequential: %f : #Trials=%ld : Elapsed-time-total(ms)=%.2f\n", pi_seqential, (long) NUM_TRIALS, comp_time_total);
}

void execute_sequential_version(void)
{
    int 	seq_id = 1;
    long	i, sequential_hits = 0;
	double	x, y;
	pi_seqential = 0;
	for (i=0; i < NUM_TRIALS; i++) {
		x = ((double) rand_r(&seq_id)) / RAND_MAX;
		y = ((double) rand_r(&seq_id)) / RAND_MAX;
		if (x*x + y*y < 1.0)
			sequential_hits++;
	}
	pi_seqential = 4.0 * sequential_hits / NUM_TRIALS;
}

//-----------------parallel code related------------------

void handle_parallel_version()
{
    struct timespec 	start_time, end_time;
    float 			    comp_time_total; 	// in milli seconds
    unsigned long 		sec, nsec;

	GET_TIME(start_time);
	execute_parallel_version();
    GET_TIME(end_time);

    comp_time_total = elapsed_time_msec(&start_time, &end_time, &sec, &nsec);

    printf("Pi-Parallel: %f : #Trials=%ld : #Threads=%ld : Elapsed-time-total(ms)=%.2f\n", pi_parallel, (long) NUM_TRIALS, num_threads, comp_time_total);
}

void execute_parallel_version()
{
    pthread_t 		    tid[MAX_THREADS];
    int 			    i, myid[MAX_THREADS];
    pi_parallel = 0;

    pthread_mutex_init(&mutex, NULL);
    for (i = 0; i < num_threads; i++) {
	    myid[i] = i;
		pthread_create(&tid[i], NULL, run, &myid[i]);
    }
    for (i = 0; i < num_threads; i++)  {
	    pthread_join(tid[i], NULL);
    }
    pi_parallel = 4.0 * parallel_hits / NUM_TRIALS;
	pthread_mutex_destroy(&mutex);
}

void *run(void * tid)
{
    int 	myid = *((int *)tid);
    long 	start = (myid * (long)NUM_TRIALS)/num_threads;
    long 	end = ((myid+1) * (long)NUM_TRIALS) /num_threads;
    long	i, my_hits = 0;
	double	x, y;

	for (i = start; i < end; i++) {
		x = ((double) rand_r(&myid)) / RAND_MAX;
		y = ((double) rand_r(&myid)) / RAND_MAX;
		if (x*x + y*y < 1.0)
	        my_hits++;
	}

   	pthread_mutex_lock (&mutex);
   	parallel_hits += my_hits;
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