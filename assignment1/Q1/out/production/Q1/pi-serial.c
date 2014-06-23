#include <stdio.h>
#include <stdlib.h>
#include <time.h>		// for clock_gettime()
#include <errno.h>		// for perror()

#define NUM_TRIALS	(100*1000*1000)
#define GET_TIME(x);	if (clock_gettime(CLOCK_MONOTONIC, &(x)) < 0) \
				{ perror("clock_gettime( ):"); exit(EXIT_FAILURE); }

long	num_hits=0;

void	hit(void);
float 	elapsed_time_msec(struct timespec *, struct timespec *, long *, long *);

int main (int argc, char *argv[])
{
	float 			pi;
	struct timespec 	t0, t1;
	unsigned long 		sec, nsec;
	float 			comp_time; 	// in milli seconds
	
	srandom(0);
	GET_TIME(t0);
	hit( );
	pi = 4.0 * num_hits / NUM_TRIALS;
	GET_TIME(t1);
	comp_time = elapsed_time_msec(&t0, &t1, &sec, &nsec);
	printf("Pi-serial: Pi = %f : #Trials=%ld : Elapsed-time(ms)=%.2f\n", pi, (long) NUM_TRIALS, comp_time);
	return 0;
}

void hit(void)
{
	long	i;
	double	x, y;
	for (i=0; i < NUM_TRIALS; i++) {
		x = ((double) random())/ RAND_MAX;
		y = ((double) random())/ RAND_MAX;
		if (x*x + y*y < 1.0) 
			num_hits++;
	}
}

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
