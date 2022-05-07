#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "helper.h"

// declarations
int substitution_matrix[][4] = {{2, 1, -1, -1}, {1, 2, -1, -1}, {-1, -1, 2, 1}, {-1, -1, 1, 2}};
void NeedlemanWunsch(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dna, char *dna_align, int tileSz);
void NeedlemanWunschTiled(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dnaseq1, char *dnaseq2, int tileSz);

// DNA substitution matrix
//    C  T  A  G
// C  2  1 -1 -1
// T  1  2 -1 -1
// A -1 -1  2  1
// G -1 -1  1  2

int main( int argc, char *argv[] )
{
    struct timeval calc;
    double calctime;
    srand((unsigned)time(NULL));

    char *dnaseq1, *dnaseq2;
    int seq_len = 5000, tileSz = 1;
    if ( argc > 1)
    {
        seq_len = atoi(argv[1]);
    }

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
    NeedlemanWunschTiled(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2, tileSz);
    calctime = tock(&calc);
    
    printf("%s\n%s\n", dnaseq1, dnaseq2);
    // for (int i = 0; i < seq_len; i++)
    // {
        // for (int j = 0; j < seq_len; j++)
        // {
            // printf("%7d ", scoringMatrix[i][j]);
        // }
        // printf("\n");
    // }

    double mem_bw = (double)seq_len * (double)seq_len * (double)sizeof(int);
    double gflops = (double)seq_len * (double)seq_len * (double)2 / calctime;

    printf("Time (in milli-secs) %f\n", calctime * 1000);
    printf("Memory Bandwidth (in GBytes/s): %f\n", mem_bw);
    printf("Compute Throughput (in GFlops/s): %f\n", gflops / calctime);

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
void NeedlemanWunsch(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dnaseq1, char *dnaseq2, int tileSz)
{
    //tiling 
    for( int sum = 0; sum < 2*n - tileSz; sum += tileSz )
    {
        int l_bound = max2(0, sum - n);
        int u_bound = min(sum, n);
        for( int ti = l_bound; ti < u_bound; ti += tileSz )
        {
            int tj = sum - ti;
            // ti, tj are the starting indices of the tileSz * tileSz tile
            // anti-diagonals
            for( int s2 = sum; s2 < sum + 2*tileSz - 1; s2++ )
            {
                int lower_bound = max2(tj, s2 - tileSz);
                int upper_bound = min(tj + tileSz, s2);
                for( int j = lower_bound; j < upper_bound; j++ )
                {
                    // if( j == 0 || s2 - j == 0)
                        // continue;
                    printf("first i, j: %d %d\n", j, s2 - j);
                    // int i = s2 - j;
                    scoringMatrix[s2 - j][j] = max(scoringMatrix[s2 - j - 1][j - 1] + substitution_matrix[ ctag(dnaseq1[s2 - j]) ][ ctag(dnaseq2[j]) ], // mismatch
                                              scoringMatrix[s2 - j - 1][j] + gap_penalty, scoringMatrix[s2 - j][j - 1] + gap_penalty);   // gap penalties
                }
            }
        }
    }

    //Tiling with anti-diagonals for each tile

}

void NeedlemanWunschTiled(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dnaseq1, char *dnaseq2, int tileSz)
{
    int total = (n / tileSz) * 2 + 1;
    for (int count = 0; count < total; count++)
    {
        int i, j;
        if (count <= (n / tileSz))
        {
            i = 0;
            j = count * tileSz;
        }
        else
        {
            i = (count - (n / tileSz)) * tileSz;
            j = (n / tileSz) * tileSz;
        }
        while (j >= 0 && i <= (n / tileSz) * tileSz)
        {
            // compute for each tile here
            // printf("i: %d, j: %d\n", i, j);
            int rowupperlimit,colupperlimit;
            rowupperlimit = min(i + tileSz, n);
            colupperlimit = min(j + tileSz, n);
            // printf("i: %d, j: %d, rowupperlimit: %d, colupperlimit: %d\n", i, j, rowupperlimit, colupperlimit);
            for(int sum = 0 ; sum < rowupperlimit-i+colupperlimit-j+1 ; sum++)
            {
                int k1,k2;
                if(sum < colupperlimit-j+1)
                {
                    k1 = i;
                    k2 = j + sum;
                }
                else
                {
                    k1 = i + sum - (colupperlimit-j);
                    k2 = colupperlimit;
                }

                int u_bound = min(i + sum, rowupperlimit);

                #pragma omp parallel for
                for( ;k1 < u_bound; k1++ )
                {
                    // printf("k1: %d, k2: %d\n", k1, k2);
                    k2 = (i + j + sum) - k1;
                    if (k1 == 0 && k2 == 0)
                    {
                        scoringMatrix[k1][k2] = substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])];
                    }
                    else if (k1 > 0 && k2 > 0)
                    {
                        scoringMatrix[k1][k2] = max(scoringMatrix[k1 - 1][k2 - 1] + substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])],
                                                    scoringMatrix[k1 - 1][k2] + gap_penalty, scoringMatrix[k1][k2 - 1] + gap_penalty);
                    }
                }
            }

            j = j - tileSz;
            i = i + tileSz;
        }
    }
}