#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "helper.h"

// declarations
int substitution_matrix[][4] = {{2, -1, -1, -1}, {-1, 2, -1, -1}, {-1, -1, 2, -1}, {-1, -1, -1, 2}};
void NeedlemanWunsch(int **scoringMatrix, int n, char *bases, int gap_penality, char *dna, char *dna_align);
void backtrack(int **scoringMatrix, int n, char *bases, int gap_penality, char *dnaseq1, char *dnaseq2);
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
    int seq_len = 5000;
    if ( argc > 1)
    {
        seq_len = atoi(argv[1]);
    }
    int tileSz = 256;
    
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
    scoringMatrix[0][0] = substitution_matrix[ctag(dnaseq1[0])][ctag(dnaseq2[0])];
    // tick(&calc);

    for (int i = 1; i < seq_len; i++)
    {
       // printf("%d ", scoringMatrix[0][i]);
        scoringMatrix[0][i] = max(scoringMatrix[0][i - 1] + gap_penalty, (i*gap_penalty) + substitution_matrix[ctag(dnaseq1[0])][ctag(dnaseq2[i])], ((i+1)*gap_penalty) + gap_penalty);
        scoringMatrix[i][0] = max(scoringMatrix[i - 1][0] + gap_penalty, (i*gap_penalty) + substitution_matrix[ctag(dnaseq1[i])][ctag(dnaseq2[0])], ((i+1)*gap_penalty) + gap_penalty);
    }
    
    for (int i = 1; i < seq_len; i++)
    {
       // printf("%d ", scoringMatrix[0][i]);
        scoringMatrix[0][i] = max(scoringMatrix[0][i - 1] + gap_penalty, (i*gap_penalty) + substitution_matrix[ctag(dnaseq1[0])][ctag(dnaseq2[i])], ((i+1)*gap_penalty) + gap_penalty);
        scoringMatrix[i][0] = max(scoringMatrix[i - 1][0] + gap_penalty, (i*gap_penalty) + substitution_matrix[ctag(dnaseq1[i])][ctag(dnaseq2[0])], ((i+1)*gap_penalty) + gap_penalty);
    }

    tick(&calc);
    NeedlemanWunschTiled(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2, tileSz);
    calctime = tock(&calc);

    backtrack(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2);

    double gflops = 3 * (double)seq_len * (double)seq_len / calctime / 1e9;
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
    for(int count = 0; count < 2*n-1; count++)
    {
        int k1;
        if(count < n)
        {
            k1 = 0;
        }
        else
        {
            k1 = count - n + 1;
        }
        // k2 >= 0 => k1<=count , k1<=n-1
        for(;k1<=min(count,n-1);k1++)
            {
                int k2 = count - k1;
                if (k1 == 0 && k2 == 0)
                {
                    scoringMatrix[k1][k2] = substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])];
                }
                else if(k1>0 && k2>0)
                {
                    scoringMatrix[k1][k2] = max(scoringMatrix[k1-1][k2] + gap_penalty, scoringMatrix[k1][k2-1] + gap_penalty, scoringMatrix[k1-1][k2-1] + substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])]);
                }
            }
    }
}
void reverse(char *s)
{
    int i, j;
    char c;
    for (i = 0, j = strlen(s) - 1; i < j; i++, j--)
    {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}
// backtracking in needleman algorithm
void backtrack(int **scoringMatrix, int n, char *bases, int gap_penality, char *dnaseq1, char *dnaseq2)
{
    int i = n - 1;
    int j = n - 1;
    char *dnaseq1_align = (char *)malloc(100000 * sizeof(char));
    char *dnaseq2_align = (char *)malloc(100000 * sizeof(char));
    int k = 0;
    while (i >= 0 && j >= 0)
    {
        // printf("i: %d, j: %d\n", i, j);
        if (i == 0 && j == 0)
        {
            dnaseq1_align[k] = dnaseq1[i];
            dnaseq2_align[k] = dnaseq2[j];
            k++;
            i--;
            j--;
        }
        // printf("%d %d %d\n",scoringMatrix[i][j],scoringMatrix[i-1][j-1],substitution_matrix[ ctag(dnaseq1[i]) ][ ctag(dnaseq2[j]) ]);
        else if (i > 0 && j > 0 && scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + substitution_matrix[ctag(dnaseq1[i])][ctag(dnaseq2[j])])
        {
            dnaseq1_align[k] = dnaseq1[i];
            dnaseq2_align[k] = dnaseq2[j];
            i--;
            j--;
            k++;
        }
        else if (i > 0)
        {
            dnaseq1_align[k] = dnaseq1[i];
            dnaseq2_align[k] = '-';
            i--;
            k++;
        }
        else if (j > 0)
        {
            dnaseq1_align[k] = '-';
            dnaseq2_align[k] = dnaseq2[j];
            j--;
            k++;
        }
    }
    while (i >= 0)
    {
        dnaseq1_align[k] = dnaseq1[i];
        dnaseq2_align[k] = '-';
        i--;
        k++;
    }
    while (j >= 0)
    {
        dnaseq2_align[k] = dnaseq2[j];
        dnaseq1_align[k] = '-';
        j--;
        k++;
    }
    reverse(dnaseq1_align);
    reverse(dnaseq2_align);
    // printf("%s %s\n", dnaseq1_align, dnaseq2_align);
}


void NeedlemanWunschTiled(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dnaseq1, char *dnaseq2, int tileSz)
{
    int total = (n / tileSz) * 2 + 1;
    for (int count = 0; count < total; count++)
    {
        int t2;
        if (count <= (n / tileSz))
        {
            t2 = 0;
          //  j = count * tileSz;
        }
        else
        {
            t2 = (count - (n / tileSz)) * tileSz;
          //  j = (n / tileSz) * tileSz;
        }
        // j >= 0 => count*tileSz - i >= 0 =>  i <= count*tileSz && i<=
        #pragma omp parallel for
        for( int i = t2; i <= min((n / tileSz) * tileSz,count*tileSz);i+=tileSz)
        {
            int j = count * tileSz - i;

            int rowupperlimit,colupperlimit;
            rowupperlimit = min(i + tileSz, n);
            colupperlimit = min(j + tileSz, n);
            // printf("i: %d, j: %d, rowupperlimit: %d, colupperlimit: %d\n", i, j, rowupperlimit, colupperlimit);


            for(int sum = 0 ; sum < rowupperlimit-i+colupperlimit-j+1 ; sum++)
            {
                int tmp1;
                if(sum < colupperlimit-j+1)
                {
                    tmp1 = i;
                }
                else
                {
                    tmp1 = i + sum - (colupperlimit-j);
                }
                // k1+k2 = i+j+sum k1 < rowupperlimit i+j+sum - k1 >= j  , k1 <= i+sum

                #pragma omp parallel for
                for( int k1 = tmp1;k1 <= min(i+sum,rowupperlimit-1);k1++)
                {
                    int k2 = i+j+sum-k1;
                    // printf("k1: %d, k2: %d\n", k1, k2);
                    if (k1 == 0 && k2 == 0)
                    {
                       // printf("k1: %d, k2: %d\n", k1, k2);
                        scoringMatrix[k1][k2] = substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])];
                    }
                    else if (k1 > 0 && k2 > 0)
                    {
                        //printf("k1: %d, k2: %d\n", k1, k2);
                        scoringMatrix[k1][k2] = max(scoringMatrix[k1 - 1][k2 - 1] + substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])],
                                                    scoringMatrix[k1 - 1][k2] + gap_penalty, scoringMatrix[k1][k2 - 1] + gap_penalty);
                    }
                }
            }
        }
    }
}