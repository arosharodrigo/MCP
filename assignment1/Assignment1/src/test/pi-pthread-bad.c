// pi-pthread-bad.c
// by Sanath Jayasena, June 2014
// Provided as material for Assignment 1
// CS5270 MCP (MSc in CS, UoM)
//
// Note: This implementation has a performance bug. The performance actually
// decreases when the number of threads increase. Provided only for demo purpose. 
//
// compile as: gcc [-DNUM_TRIALS=x] -lrt -lpthread

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>		// for clock_gettime()
#include <errno.h>		// for perror()

#define MAX_THREADS     32
#define NUM_TRIALS	(100*1000*1000)
#define GET_TIME(x);	if (clock_gettime(CLOCK_MONOTONIC, &(x)) < 0) \
				{ perror("clock_gettime( ):"); exit(EXIT_FAILURE); }

void	*hit(void *);
float 	elapsed_time_msec(struct timespec *, struct timespec *, long *, long *);

long    num_hits=0, num_threads=0;
pthread_mutex_t mutex;


int main (int argc, char *argv[])
{
	float 			pi;
	pthread_t 		tid[MAX_THREADS];
	int 			i, *myid[MAX_THREADS];
	struct timespec 	t0, t1;
	unsigned long 		sec, nsec;
	float 			comp_time; 	// in milli seconds
	
    	if (argc != 2) {
		printf("Usage: %s <num_threads>\n", argv[0]);
		exit(0);
    	}
    	num_threads = atoi(argv[1]);
    	if (num_threads > MAX_THREADS) {
		printf("num_threads > MAXTHREADS (%d)\n", MAX_THREADS);
		exit(0);
    	}
	GET_TIME(t0);
	pthread_mutex_init(&mutex, NULL);
    	for (i = 0; i < num_threads; i++) {                  
		*myid[i] = i;
		pthread_create(&tid[i], NULL, hit, (void *) myid[i]);
    	}                                                 
    	for (i = 0; i < num_threads; i++)  {
		pthread_join(tid[i], NULL);
    	}
	pthread_mutex_destroy(&mutex);
	pi = 4.0 * num_hits / NUM_TRIALS;
	GET_TIME(t1);
	comp_time = elapsed_time_msec(&t0, &t1, &sec, &nsec);
	printf("Pi-parallel: #Threads=%ld : Pi = %f : #Trials=%ld : Elapsed-time(ms)=%.2f\n", \
		num_threads, pi, (long) NUM_TRIALS, comp_time);
	/* Last thing that main() should do */
	pthread_exit(NULL);
}
//-----------------------------------------------
void *hit(void * tid)
{
    int 	myid = *((int *)tid);
    long 	start = (myid * (long)NUM_TRIALS)/num_threads;
    long 	end = ((myid+1) * (long)NUM_TRIALS) /num_threads;
    long	i, my_hits=0;
    struct  timespec	t10, t11;
    float   comp_time; 	// in milli seconds
    unsigned long 		sec, nsec;
    double	x, y;

    GET_TIME(t10);
	printf("Thread %d: started at %.2f\n", myid, &t10);

	for (i=start; i < end; i++) {
		x = ((double) random())/ RAND_MAX;
		y = ((double) random())/ RAND_MAX;
		if (x*x + y*y < 1.0)
			my_hits++;
	}

//   	pthread_mutex_lock (&mutex);
   	num_hits += my_hits;
//   	pthread_mutex_unlock (&mutex);
   	GET_TIME(t11);
   	comp_time = elapsed_time_msec(&t10, &t11, &sec, &nsec);
	printf("Thread %d: start=%ld, end=%ld, hits=%ld, Elapsed-time(ms)=%.2f\n", myid, start, end, my_hits, comp_time);
   	pthread_exit(NULL);	
}
//-----------------------------------------------
float elapsed_time_msec(struct timespec *begin, struct timespec *end, long *sec, long *nsec)
{
	if (end->tv_nsec < begin->tv_nsec) {
		*nsec = 1000000000 - (begin->tv_nsec - end->tv_nsec);
    	*sec  = end->tv_sec - begin->tv_sec -1;
    }
    else {
	    *nsec = end->tv_nsec - begin->tv_nsec;
    	*sec  = end->tv_sec - begin->tv_sec;
	}
	return (float) (*sec) * 1000 + ((float) (*nsec)) / 1000000;

}
