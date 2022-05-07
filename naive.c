#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "helper.h"

#define DEFAULT_LEN 5000
// declarations
int substitution_matrix[][4] = {{2, 1, -1, -1}, {1, 2, -1, -1}, {-1, -1, 2, 1}, {-1, -1, 1, 2}};
void NeedlemanWunsch(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dna, char *dna_align);

// DNA substitution matrix
//    C  T  A  G
// C  2  1 -1 -1
// T  1  2 -1 -1
// A -1 -1  2  1
// G -1 -1  1  2

int main(int argc, char *argv[])
{
    struct timeval calc;
    double calctime;
    srand((unsigned)time(NULL));

    char *dnaseq1, *dnaseq2;
    int seq_len = DEFAULT_LEN;
    if (argc > 1)
        seq_len = atoi(argv[1]);

    dnaseq1 = (char *)malloc(seq_len * sizeof(char));
    dnaseq2 = (char *)malloc(seq_len * sizeof(char));
    char bases[] = {'C', 'T', 'A', 'G'};
    RandomString(seq_len, dnaseq1, bases);
    RandomString(seq_len, dnaseq2, bases);

    int gap_penalty = -1;

    int **scoringMatrix = (int **)malloc(seq_len * sizeof(int *));
    for (int i = 0; i < seq_len; i++)
    {
        scoringMatrix[i] = (int *)malloc(seq_len * sizeof(int));
    }

    for (int i = 0; i < seq_len; i++)
    {
        scoringMatrix[i][0] = substitution_matrix[ctag(bases[0])][ctag(bases[i])];
        scoringMatrix[0][i] = substitution_matrix[ctag(bases[i])][ctag(bases[0])];
    }

    tick(&calc);
    NeedlemanWunsch(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2);
    calctime = tock(&calc);
    
    // printf("%s\n%s\n", dnaseq1, dnaseq2);
    // for (int i = 0; i < seq_len; i++)
    // {
        // for (int j = 0; j < seq_len; j++)
        // {
            // printf("%7d ", scoringMatrix[i][j]);
        // }
        // printf("\n");
    // }

    // double mem_bw = (double)seq_len * (double)seq_len * (double)sizeof(int);
    double gflops = (3 * seq_len * seq_len / calctime) * 1.0e-9;

    // printf("%f\n", calctime * 1000);
    // printf("Memory Bandwidth (in GBytes/s): %f\n", mem_bw);
    printf("%lf\n", gflops);

    // free everything
    free(dnaseq1);
    free(dnaseq2);
    for (int i = 0; i < seq_len; i++)
    {
        free(scoringMatrix[i]);
    }
    free(scoringMatrix);

    return 0;
}

// Needleman-Wunsch algorithm for global alignment
void NeedlemanWunsch(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dnaseq1, char *dnaseq2)
{
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
        {
            // printf("j: %d, k: %d\n",j, i-j);
            scoringMatrix[i][j] = max(scoringMatrix[i - 1][j - 1] + substitution_matrix[ ctag(dnaseq1[i]) ][ ctag(dnaseq2[j]) ], // mismatch
                                          scoringMatrix[i - 1][j] + gap_penalty, scoringMatrix[i][j - 1] + gap_penalty);   // gap penalties
        }
    }
}