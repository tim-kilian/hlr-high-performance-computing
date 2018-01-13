/*
This by now is just a roadmap as to how it should be done
*/

#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "partdiff-par.h"
#include <sys/time.h>

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

struct timeval start_time;
struct timeval comp_time;

static void initVariables(struct calculation_arguments *arguments, struct calculation_results *results, struct options *options) {
    arguments->N = options->interlines * 8 + 9 - 1;
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
    arguments->h = (float) (1.0 / arguments->N);

    results->m = 0;
    results->stat_iteration = 0;
    results->stat_precision = 0;
}


int
calculate_lines(int myrank, int size, int interlines)
{
    int lines;
    int mylines;

    lines = 8*interlines + 9 - 2;
    
    if ((lines % size) == 0){
        mylines = lines/size;
    }   
    else {
        /* Der Rest wird auf die letzten Prozesse aufgeteilt */
        if (myrank >= (size - (lines % size)) ){
            mylines = lines/size + 1;
        }
        else {
            mylines = lines/size;
        }
    }

    return mylines;
}



static void *allocateMemory(size_t size) {
    void *p;

    if ((p = malloc(size)) == NULL) {
        printf("\n\nSpeicherprobleme!\n");
        exit(1);
    }

    return p;
}



static void allocateMatrices(int myrank, int size, struct calculation_arguments *arguments, struct options *options) {
    int N = arguments->N;
    int interlines = options->interlines;
    int lines;
    lines = calculate_lines(myrank, size, interlines)+2;
    arguments->M = allocateMemory(arguments->num_matrices * lines * (N + 1) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double **));
    for (int i = 0; i < arguments->num_matrices; i++) {
        arguments->Matrix[i] = allocateMemory(lines * sizeof(double *));
        for (int j = 0; j < lines; j++) {
            arguments->Matrix[i][j] = (double *) (arguments->M + (i * lines * (N + 1)) + (j * (N + 1)));
        }
    }
}



/* Matrix initialisieren - Alloziieren muss schon geschehen sein */
static void initMatrices(int myrank , int size, struct calculation_arguments *arguments,  struct options *options) {
    double ***Matrix = arguments->Matrix;
    int N = arguments->N;
    double h = arguments->h;
    int lines, start ;
    int interlines = options->interlines;
    lines = calculate_lines(myrank, size, interlines)+2;
    /* initialize matrix/matrices with zeros */
    for (int g = 0; g < arguments->num_matrices; g++) {
        for (int i = 0; i < lines ; i++) {
            for (int j = 0; j <= N; j++) {
                Matrix[g][i][j] = 0;
            }
        }
    }
    if (options->inf_func == FUNC_F0) {
        for (int g = 0; g < arguments->num_matrices; g++) {
            for (int j = 0; j <= N; j++){
            if (myrank == 0){
            Matrix[g][0][j] = 1 - (h * j);
            }
        if (myrank == size -1){
            Matrix[g][lines-1][j]= (h * j);
        }
            }    
        }
    }

    start = 0;
    for (int i = 0; i < myrank; i++){
        start = start + calculate_lines(i,size, options->interlines);
    }
   

    if (options->inf_func == FUNC_F0){
        for (int g = 0; g < arguments->num_matrices; g++) {
             for (int t = 0; t < lines ; t++){
         Matrix[g][t][0] = 1 - ((h * t)+(start * h));
         Matrix[g][t][N] = ((h * t) + (start*h));   
             }
    }
    }
}


/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (int myrank, int size, struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
        int i, j, test;                                   /* local variables for loops */
        int m1, m2;                                 /* used as indices for old and new matrices */
/*      int numthreads, threadnum, low, high; */
        double star;                                /* four times center value minus 4 neigh.b values */
        double residuum;                            /* residuum of current iteration */
        double maxresiduum;                         /* maximum residuum value of a slave in iteration */
    int lines;
        int const N = arguments->N;
        double const h = arguments->h;
        int start;
        double pih = 0.0;
        double fpisin = 0.0;

        int term_iteration = options->term_iteration;
        lines = calculate_lines(myrank,size, options->interlines);
    test=99;
        /* initialize m1 and m2 depending on algorithm */
        if (options->method == METH_JACOBI)
        {
                m1 = 0;
                m2 = 1;
        }
        else
        {
                m1 = 0;
                m2 = 0;
        }


        if (options->inf_func == FUNC_FPISIN)
        {
                pih = PI * h;
                fpisin = 0.25 * TWO_PI_SQUARE * h * h;
        }

    start = 0;
        for (int i = 0; i < myrank; i++){
            start = start + calculate_lines(i,size, options->interlines);
        }
        while (term_iteration > 0)
        {
                double** Matrix_Out = arguments->Matrix[m1];
                double** Matrix_In  = arguments->Matrix[m2];

                maxresiduum = 0;

        if (myrank < (size-1)){
                MPI_Send(Matrix_In[lines], N+1 , MPI_DOUBLE, myrank + 1, 99, MPI_COMM_WORLD);
                MPI_Recv(Matrix_In[lines + 1], N+1, MPI_DOUBLE, myrank + 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (myrank > 0){
                MPI_Recv(Matrix_In[0], N + 1, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(Matrix_In[1], N+1, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD);
        }


                for (i = 1; i <= lines; i++)
                {
                        double fpisin_i = 0.0;
                    double position = i + (double) start;
                        if (options->inf_func == FUNC_FPISIN)
                        {
                                fpisin_i = fpisin * sin(pih * (double)position);
                        }

                        /* over all columns */

                        for (j = 1; j < N; j++)
                        {
                                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

                                if (options->inf_func == FUNC_FPISIN)
                                {
                                        star += fpisin_i * sin(pih * (double)j);
                                }
                                if (options->termination == TERM_PREC || term_iteration == 1)
                                {
                                        residuum = Matrix_In[i][j] - star;
                                        residuum = (residuum < 0) ? -residuum : residuum;
                                        maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                                }

                                Matrix_Out[i][j] = star;
                        }
                }
        

                results->stat_iteration++;
                results->stat_precision = maxresiduum;



                /* exchange m1 and m2 */
                i = m1;
                m1 = m2;
                m2 = i;



                /* check for stopping calculation depending on termination method */
                if (options->termination == TERM_PREC)
                {
                        if (maxresiduum < options->term_precision)
                        {
                                term_iteration = 0;
                        }
                MPI_Allreduce(&term_iteration, &test, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                term_iteration = test;
                }
                else if (options->termination == TERM_ITER)
                {
                        term_iteration--;
                }
        }
        results->m = m2;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_gs (int myrank, int size, struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
        int i, j, test;                                   /* local variables for loops */
        int m1, m2;                                 /* used as indices for old and new matrices */
/*      int numthreads, threadnum, low, high; */
        double star;                                /* four times center value minus 4 neigh.b values */
        double residuum;                            /* residuum of current iteration */
        double maxresiduum, buffer_maxresiduum;                         /* maximum residuum value of a slave in iteration */
    int lines;
        int const N = arguments->N;
        double const h = arguments->h;
        int start;
    int iteration;
        double pih = 0.0;
        double fpisin = 0.0;
    MPI_Status status;
        iteration = 0;
        int a = 0;

    int help;
    help = 0;

        int term_iteration = options->term_iteration;
        lines = calculate_lines(myrank,size, options->interlines);
    test=99;
        /* initialize m1 and m2 depending on algorithm */
        if (options->method == METH_JACOBI)
        {
                m1 = 0;
                m2 = 1;
        }
        else
        {
                m1 = 0;
                m2 = 0;
        }


        if (options->inf_func == FUNC_FPISIN)
        {
                pih = PI * h;
                fpisin = 0.25 * TWO_PI_SQUARE * h * h;
        }

    start = 0;
        for (int i = 0; i < myrank; i++){
            start = start + calculate_lines(i,size, options->interlines);
        }
    term_iteration = term_iteration + myrank;
        while (term_iteration > 0)
        {
                double** Matrix_Out = arguments->Matrix[m1];
                double** Matrix_In  = arguments->Matrix[m2];
        if ((options->termination == TERM_PREC) &&(help == 1)){
            term_iteration--;
        }
                maxresiduum = 0;

        if (myrank < (size-1)){
            MPI_Send(Matrix_In[lines], N+1, MPI_DOUBLE, myrank +1, myrank, MPI_COMM_WORLD);
        }
        if (myrank > 0){
                        if (term_iteration > 1){
                MPI_Recv(Matrix_In[0], N+1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD, &status);
            }
        }

        
        iteration++;
        if (myrank < iteration){

                for (i = 1; i <= lines; i++)
                {
                        double fpisin_i = 0.0;
                    double position = i + (double) start;
                        if (options->inf_func == FUNC_FPISIN)
                        {
                                fpisin_i = fpisin * sin(pih * (double)position);
                        }

                        /* over all columns */

                        for (j = 1; j < N; j++)
                        {
                                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

                                if (options->inf_func == FUNC_FPISIN)
                                {
                                        star += fpisin_i * sin(pih * (double)j);
                                }
                                if (options->termination == TERM_PREC || term_iteration == 1)
                                {
                                        residuum = Matrix_In[i][j] - star;
                                        residuum = (residuum < 0) ? -residuum : residuum;
                                        maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                                }

                                Matrix_Out[i][j] = star;
                        }

      if ((myrank > 0) && (i==1))
        {
           if (term_iteration > 1)
            {
            MPI_Send(Matrix_In[1], N+1, MPI_DOUBLE, myrank -1, myrank, MPI_COMM_WORLD);
            }
        }

        if ( (myrank < (size-1)) && (i==1) && (myrank+1 < iteration)){
            MPI_Recv(Matrix_In[lines+1], N+1, MPI_DOUBLE, myrank+1, myrank +1, MPI_COMM_WORLD, &status);
            a += 1;
        }
        
                }
        }


                results->stat_iteration++;
                results->stat_precision = maxresiduum;



                /* exchange m1 and m2 */
                i = m1;
                m1 = m2;
                m2 = i;

        
                /* check for stopping calculation depending on termination method */
                if ((options->termination == TERM_PREC) && (iteration > size) && (help == 0))
                {
                        if (maxresiduum < options->term_precision)
                        {
                                term_iteration = 0;
                        }
                MPI_Allreduce(&term_iteration, &test, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allreduce(&maxresiduum, &buffer_maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                maxresiduum = buffer_maxresiduum;
                term_iteration = test; // ist 0 wenn alle Praezision erreicht haben
            if (term_iteration == 0){
                help = 1;
                term_iteration = term_iteration + myrank;
            }
                }
                else if (options->termination == TERM_ITER)
                {
                        term_iteration--;
                }
        // In der letzten Iteration muss ebenfalls gesendet werden, deshalb hier erneut MPI_Send/Recv
        if ((term_iteration == 0) && (myrank < (size -1))){
            MPI_Send(Matrix_In[lines], N+1, MPI_DOUBLE, myrank +1, myrank, MPI_COMM_WORLD);
        }
        if ((term_iteration == 1) && (myrank > 0)){
            MPI_Recv(Matrix_In[0], N+1, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD, &status);
        }

        }
        results->m = m2;
}


/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics(struct calculation_arguments const *arguments, struct calculation_results const *results, struct options const *options) {
    int N = arguments->N;
    double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

    printf("Berechnungszeit:    %f s \n", time);
    printf("Speicherbedarf:     %f MiB\n",
           (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
    printf("Berechnungsmethode: ");

    if (options->method == METH_GAUSS_SEIDEL) {
        printf("Gauß-Seidel");
    } else if (options->method == METH_JACOBI) {
        printf("Jacobi");
    }

    printf("\n");
    printf("Interlines:         %d\n", options->interlines);
    printf("Stoerfunktion:      ");

    if (options->inf_func == FUNC_F0) {
        printf("f(x,y) = 0");
    } else if (options->inf_func == FUNC_FPISIN) {
        printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
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
    printf("\n");
}


static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;
  double** Matrix = arguments->Matrix[results->m];
  MPI_Status status;
  MPI_Barrier(MPI_COMM_WORLD);
  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;

  /* last line belongs to rank size - 1 */
  if (rank + 1 == size)
    to++;

  if (rank == 0)
    printf("Matrix:\n");

  for (y = 0; y < 9; y++)
  {
    int line = y * (options->interlines + 1);

    if (rank == 0)
    {
      //printf("iteration y: %d\n", y);
      /* check whether this line belongs to rank 0 */
      if (line < from || line > to)
      {
        /* use the tag to receive the lines in the correct order
         * the line is stored in Matrix[0], because we do not need it anymore */
        MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }

    if (rank == 0)
    {
      for (x = 0; x < 9; x++)
      {
        int col = x * (options->interlines + 1);
        if (line >= from && line <= to)
        {
          /* this line belongs to rank 0 */
          printf("%7.4f", Matrix[line][col]);
        }
        else
        {
          /* this line belongs to another rank and was received above */
          printf("%7.4f", Matrix[0][col]);
        }
      }

      printf("\n");
    }
  }

  fflush(stdout);
}



int 
main(int argc, char** argv)
{
    MPI_Init( &argc, &argv );
    int myrank, size, from, to, cancel;
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;
    cancel = 0;
    from = 0;
    to = 0;
    // get parameters
    AskParams(&options, argc, argv);
    initVariables(&arguments, &results, &options);
    /* Parallel starts here */
    /* MPI Initialisierung, Rangabfrage und Anzahlabfrage */
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    //  get and initialize variables and matrices
    allocateMatrices(myrank, size, &arguments, &options);
    initMatrices(myrank , size, &arguments, &options);

    gettimeofday(&start_time, NULL);

    if (options.method == METH_JACOBI){
        calculate(myrank, size, &arguments, &results, &options);
    }
    else {
        for (int i = 0; i < size; i++){
                cancel = calculate_lines(i,size, options.interlines);
            if (cancel == 1){
                if (myrank == 0) {printf("For this Gauss-Seidel implementation each Process needs at least 2 lines for calculations.\n");}
                exit(0);
            }
            }

        calculate_gs(myrank, size, &arguments, &results, &options);
    }

    gettimeofday(&comp_time, NULL);
    
    if (myrank == 0){
    displayStatistics(&arguments, &results, &options);
    }

    for (int i = 0; i < myrank; i++){
    from = from + calculate_lines(i,size, options.interlines);
    }
    from = from + 1;

   for (int i = 0; i <= myrank; i++){
    to = to + calculate_lines(i, size, options.interlines);
    }
  
    DisplayMatrix(&arguments, &results, &options, myrank, size, from, to);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}