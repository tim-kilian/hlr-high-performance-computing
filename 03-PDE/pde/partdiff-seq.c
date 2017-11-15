/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**            JK und andere besseres Timing, FLOP-Berechnung              **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
/**            Jacobi methods.                                             **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include "partdiff-seq.h"

struct calculation_arguments {
    int N;              /* number of spaces between lines (lines=N+1)     */
    int num_matrices;   /* number of matrices                             */
    double ***Matrix;      /* index matrix used for addressing M             */
    double *M;             /* two matrices with real values                  */
    double h;              /* length of a space between two lines            */
};

struct calculation_results {
    int m;
    int stat_iteration; /* number of current iteration                    */
    double stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static void initVariables(struct calculation_arguments *arguments, struct calculation_results *results, struct options *options) {
    arguments->N = options->interlines * 8 + 9 - 1;
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
    arguments->h = (float) (1.0 / arguments->N);

    results->m = 0;
    results->stat_iteration = 0;
    results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static void freeMatrices(struct calculation_arguments *arguments) {
    for (int i = 0; i < arguments->num_matrices; i++) {
        free(arguments->Matrix[i]);
    }

    free(arguments->Matrix);
    free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static void *allocateMemory(size_t size) {
    void *p;

    if ((p = malloc(size)) == NULL) {
        printf("\n\nSpeicherprobleme!\n");
        exit(1);
    }

    return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void allocateMatrices(struct calculation_arguments *arguments) {
    int N = arguments->N;

    arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double **));

    for (int i = 0; i < arguments->num_matrices; i++) {
        arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double *));

        for (int j = 0; j <= N; j++) {
            arguments->Matrix[i][j] = (double *) (arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1)));
        }
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices(struct calculation_arguments *arguments, struct options *options) {
    int N = arguments->N;
    double h = arguments->h;
    double ***Matrix = arguments->Matrix;

    /* initialize matrix/matrices with zeros */
    for (int g = 0; g < arguments->num_matrices; g++) {
        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= N; j++) {
                Matrix[g][i][j] = 0;
            }
        }
    }

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0) {
        for (int i = 0; i <= N; i++) {
            for (int j = 0; j < arguments->num_matrices; j++) {
                Matrix[j][i][0] = 1 - (h * i);
                Matrix[j][i][N] = h * i;
                Matrix[j][0][i] = 1 - (h * i);
                Matrix[j][N][i] = h * i;
            }
        }

        for (int j = 0; j < arguments->num_matrices; j++) {
            Matrix[j][N][0] = 0;
            Matrix[j][0][N] = 0;
        }
    }
}

/* ************************************************************************ */
/* getResiduum: calculates residuum                                         */
/* Input: x,y - actual column and row                                       */
/* ************************************************************************ */
double getResiduum(double h, struct options *options, int x, int y, double star) {
    if (options->inf_func == FUNC_F0) {
        return -star / 4.0;
    }
    return (f(x*h,y*h) * h * h - star) / 4.0;
}

double f(double xh, double yh) {
    return (double) TWO_PI_SQUARE * sin(yh * PI) * sin(xh * PI);
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void calculate(struct calculation_arguments *arguments, struct calculation_results *results, struct options *options) {
    double ***M = arguments->Matrix;

    /* initialize m1 and m2 depending on algorithm */
    int m1 = 0;
    int m2 = (options->method == METH_GAUSS_SEIDEL) ? 0 : 1;

    while (options->term_iteration > 0) {
        double max_residuum = 0;

        for (int i = 1; i < arguments->N; i++) {
            for (int j = 1; j < arguments->N; j++) {
                double star = 4.0 * M[m2][i][j] - M[m2][i - 1][j] - M[m2][i + 1][j] - M[m2][i][j - 1] - M[m2][i][j + 1];

                double korrektur = getResiduum(arguments->h, options, i, j, star);
                double residuum = fabs(korrektur);
                max_residuum = fmax(max_residuum, residuum);

                M[m1][i][j] = M[m2][i][j] + korrektur;
            }
        }

        results->stat_iteration++;
        results->stat_precision = max_residuum;

        int temp = m1; m1 = m2; m2 = temp;

        // check for stopping calculation, depending on termination method
        if (options->termination == TERM_PREC) {
            if (max_residuum < options->term_precision) options->term_iteration = 0;
        } else if (options->termination == TERM_ITER) {
            options->term_iteration--;
        }
    }

    results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics(struct calculation_arguments *arguments, struct calculation_results *results, struct options *options) {
    (void) arguments;

    double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;
    printf("Berechnungszeit:    %f s \n", time);

    printf("Berechnungsmethode: ");

    if (options->method == METH_GAUSS_SEIDEL) {
        printf("Gauss-Seidel");
    } else if (options->method == METH_JACOBI) {
        printf("Jacobi");
    }

    printf("\n");
    printf("Interlines:         %d\n", options->interlines);
    printf("Stoerfunktion:      ");

    if (options->inf_func == FUNC_F0) {
        printf("f(x,y)=0");
    } else if (options->inf_func == FUNC_FPISIN) {
        printf("f(x,y)=2pi^2*sin(pi*x)sin(pi*y)");
    }

    printf("\n");
    printf("Terminierung:       ");

    if (options->termination == TERM_PREC) {
        printf("Hinreichende Genaugkeit");
    } else if (options->termination == TERM_ITER) {
        printf("Anzahl der Iterationen");
    }

    printf("\n");
    printf("Anzahl Iterationen: %d\n", results->stat_iteration);
    printf("Norm des Fehlers:   %e\n", results->stat_precision);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int main(int argc, char **argv) {
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;

    // get parameters
    AskParams(&options, argc, argv);

    initVariables(&arguments, &results, &options);

    //  get and initialize variables and matrices
    allocateMatrices(&arguments);
    initMatrices(&arguments, &options);

    // start timer
    gettimeofday(&start_time, NULL);
    //  solve the equation
    calculate(&arguments, &results, &options);
    // stop timer
    gettimeofday(&comp_time, NULL);

    displayStatistics(&arguments, &results, &options);
    // display some
    // statistics and
    DisplayMatrix("Matrix:", arguments.Matrix[results.m][0], options.interlines);

    // free memory
    freeMatrices(&arguments);

    return 0;
}
