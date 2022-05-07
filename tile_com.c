#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// declarations
int substitution_matrix[][4] = {{2, 1, -1, -1}, {1, 2, -1, -1}, {-1, -1, 2, 1}, {-1, -1, 1, 2}};
void RandomString(int n, char *dna, char *bases);
void NeedlemanWunsch(int **scoringMatrix, int n, char *bases, int gap_penality, char *dna, char *dna_align);
void backtrack(int **scoringMatrix, int n, char *bases, int gap_penality, char *dnaseq1, char *dnaseq2);
void NeedlemanWunschTiled(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dnaseq1, char *dnaseq2, int tileSz);

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
    if (c == 'C')
    {
        return 0;
    }
    else if (c == 'T')
    {
        return 1;
    }
    else if (c == 'A')
    {
        return 2;
    }
    else
    {
        return 3;
    }
}
// DNA substitution matrix
//    C  T  A  G
// C  2  1 -1 -1
// T  1  2 -1 -1
// A -1 -1  2  1
// G -1 -1  1  2

int main()
{
    struct timeval calc;
    double calctime;
    srand((unsigned)time(NULL));

    char *dnaseq1, *dnaseq2;
    int seq_len = 8;
    dnaseq1 = (char *)malloc(seq_len * sizeof(char));
    dnaseq2 = (char *)malloc(seq_len * sizeof(char));
    char bases[] = {'C', 'T', 'A', 'G'};
    RandomString(seq_len, dnaseq1, bases);
    RandomString(seq_len, dnaseq2, bases);
    int gap_penalty = -1;
    //  int ctag[] = {'C' - 'A','T' - 'A','A' - 'A', 'G' - 'A'};

    // tick(&calc);
    // NeedlemanWunsch(dnaseq1, dnaseq2, seq_len, gap_penalty, substitution_matrix);
    // calctime = tock(&calc);

    int *scoringMatrix = (int *)malloc(seq_len * sizeof(int *));
    for (int i = 0; i < seq_len; i++)
    {
        scoringMatrix[i] = (int *)malloc(seq_len * sizeof(int));
    }

    // tick(&calc);
    for (int i = 0; i < seq_len; i++)
    {
        scoringMatrix[i][0] = substitution_matrix[ctag(bases[0])][ctag(bases[i])];
        scoringMatrix[0][i] = substitution_matrix[ctag(bases[i])][ctag(bases[0])];
    }

     NeedlemanWunsch(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2);
    //NeedlemanWunschTiled(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2, 3);
     // calctime = tock(&calc);
     printf("%s %s\n", dnaseq1, dnaseq2);
     for (int i = 0; i < seq_len; i++)
     {
         for (int j = 0; j < seq_len; j++)
         {
             printf("%d ", scoringMatrix[i][j]);
         }
         printf("\n");
     }
     backtrack(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2);
     for(int i=0;i<seq_len;i++)
     {
         for(int j=0;j<seq_len;j++)
         {
             scoringMatrix[i][j] = 0;
         }
     }
      for (int i = 0; i < seq_len; i++)
    {
        scoringMatrix[i][0] = substitution_matrix[ctag(bases[0])][ctag(bases[i])];
        scoringMatrix[0][i] = substitution_matrix[ctag(bases[i])][ctag(bases[0])];
    }

    NeedlemanWunschTiled(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2, 3);
    for (int i = 0; i < seq_len; i++)
     {
         for (int j = 0; j < seq_len; j++)
         {
             printf("%d ", scoringMatrix[i][j]);
         }
         printf("\n");
     }
     backtrack(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2);
   // NeedlemanWunschTiled(scoringMatrix, seq_len, bases, gap_penalty, dnaseq1, dnaseq2, 3);
    // printf("%s %s\n", dnaseq1, dnaseq2);
    // free everything

    return 0;
}

// Needleman-Wunsch algorithm for global alignment
void NeedlemanWunsch(int **scoringMatrix, int n, char *bases, int gap_penalty, char *dnaseq1, char *dnaseq2)
{
    for(int count = 0; count < 2*n-1; count++)
    {
        int k1,k2;
        if(count < n)
        {
            k1 = 0;
            k2 = count;
        }
        else
        {
            k1 = count - n + 1;
            k2 = n - 1;
        }
        while(k2 >= 0 && k1 < n)
            {
                if (k1 == 0 && k2 == 0)
                {
                    scoringMatrix[k1][k2] = substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])];
                }
                else if(k1>0 && k2>0)
                {
                    scoringMatrix[k1][k2] = max(scoringMatrix[k1-1][k2] + gap_penalty, scoringMatrix[k1][k2-1] + gap_penalty, scoringMatrix[k1-1][k2-1] + substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])]);
                }
                k1++;
                k2--;
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
        printf("i: %d, j: %d\n", i, j);
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
    printf("%s %s\n", dnaseq1_align, dnaseq2_align);
}

int min(int a, int b)
{
    return (a < b) ? a : b;
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
            printf("i: %d, j: %d, rowupperlimit: %d, colupperlimit: %d\n", i, j, rowupperlimit, colupperlimit);
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
                while(k2 >= j && k1 < rowupperlimit)
                {
                    // printf("k1: %d, k2: %d\n", k1, k2);
                    if (k1 == 0 && k2 == 0)
                    {
                        scoringMatrix[k1][k2] = substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])];
                    }
                    else if (k1 > 0 && k2 > 0)
                    {
                        scoringMatrix[k1][k2] = max(scoringMatrix[k1 - 1][k2 - 1] + substitution_matrix[ctag(dnaseq1[k1])][ctag(dnaseq2[k2])],
                                                    scoringMatrix[k1 - 1][k2] + gap_penalty, scoringMatrix[k1][k2 - 1] + gap_penalty);
                    }
                    k1++;
                    k2--;
                }
            }

            j = j - tileSz;
            i = i + tileSz;
        }
    }
}