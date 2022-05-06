#define _XOPEN_SOURCE

#include <sys/time.h>
#include <stdlib.h>
// #include "helper.h"

void tick(struct timeval *t)
{
    gettimeofday(t, NULL);
}

double tock(struct timeval *t)
{
    struct timeval now;
    gettimeofday(&now, NULL);
    return(double) (now.tv_sec - t->tv_sec) + 
    ((double)(now.tv_usec - t->tv_usec)/1000000.);
}

// generating random dna sequences of length n
void RandomString(int n, char *dna, char *bases)
{
    for (int i = 0; i < n; i++)
    {
        dna[i] = bases[rand() % 4];
    }
}

int max(int a, int b, int c)
{
    if (a > b)
        return (a > c ? a : c);

    return (b > c ? b : c);
}

int ctag(char c)
{
    switch (c)
    {
        case 'C':
            return 0;
        case 'T':
            return 1;
        case 'A':
            return 2;
        default:
            return 3;
    }
}