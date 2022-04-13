#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/* Standard C Function: Greatest Common Divisor */
int gcd(int a, int b)
{
	int c;
    #pragma omp shared(a, b) private(c)
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
    #pragma parallel for shared(sum, num, den) private(i, j, factor, ii, done, n)
	for (i = start; i <= end; i++){
		ii = i - start;
		sum = 1 + i;
		the_num[ii] = i;
		done = i;
		factor = 2;
        #pragma parallel reduction shared(sum)
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

int main(int argc, char **argv)
{
	long int start;
	long int end;

	scanf("%ld %ld", &start, &end);
	// printf("Number %ld to %ld\n", start, end);
	friendly_numbers(start, end);

	return EXIT_SUCCESS;
}
