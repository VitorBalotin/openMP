#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int threads = 16;

/* Standard C Function: Greatest Common Divisor */
int gcd(int a, int b)
{
	int c;
	// omp_set_num_threads(threads);
    #pragma omp shared(a, b, c)
	while (a != 0){
		c = a;
		a = b % a;
		b = c;
	}
	return b;
}

void friendly_numbers(long int start, long int end)
{
	int c = 0;
	long int last = end - start + 1;

	long int *the_num;
	the_num = (long int *)malloc(sizeof(long int) * last);
	long int *num;
	num = (long int *)malloc(sizeof(long int) * last);
	long int *den;
	den = (long int *)malloc(sizeof(long int) * last);

	long int i, j, factor, ii, sum, done, n;
	omp_set_num_threads(threads);
    #pragma omp parallel for shared(num, the_num, den, i) private(ii, done, n, sum, factor)
	for (i = start; i <= end; i++){
		ii = i - start;
		sum = 1 + i;
		the_num[ii] = i;
		done = i;
		factor = 2;
		// int nt = omp_get_num_threads();
		// printf("threads - %d\n", nt);
        #pragma omp parallel reduction(+:c) shared(sum, factor, done, i)
		while (factor < done){
			if ((i % factor) == 0){
				sum += (factor + (i / factor));
				if ((done = i / factor) == factor)
					sum -= factor;
			}
			factor++;
		}
		num[ii] = sum;
		den[ii] = i;
		n = gcd(num[ii], den[ii]);
		num[ii] /= n;
		den[ii] /= n;
	} // end for
	// omp_set_num_threads(threads);
    #pragma omp barrier
	for (i = 0; i < last; i++){
		for (j = i + 1; j < last; j++){
			if ((num[i] == num[j]) && (den[i] == den[j])){
				/* para mostrar os pares */
				// printf("%ld and %ld are friendly numbers\n", i + start, j + start);
				c++;
			}
		}
	}
	printf("Found %d pairs\n", c);

	free(the_num);
	free(num);
	free(den);
}

int main(int argc, char **argv){
	// long int start = 2;
	// long int end = 10000;
	// for(int i = 0; i < argc; i++){
	// 	printf("argumento - %d, valor é - %s.\n", i, argv[i]);
	// }
	char *p;
	long int start = strtol(argv[1], &p, 10);
	long int end = strtol(argv[2], &p, 10);
	threads = strtol(argv[3], &p, 10);

	// scanf("%ld %ld", &start, &end);
	// printf("Number %ld to %ld\n", start, end);
	friendly_numbers(start, end);

	return EXIT_SUCCESS;
}
