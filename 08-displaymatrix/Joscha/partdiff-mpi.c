/*
This by now is just a roadmap as to how it should be done
*/

#include <stdlib.h>
#include <malloc.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "partdiff-mpi.h"

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
        //printf("%d\n",g);
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
    // Fehler ab hier
    start = 0;
    for (int i = 0; i < myrank; i++){
        start = start + calculate_lines(i,size, options->interlines);
    }
   

    if (options->inf_func == FUNC_F0){
        for (int g = 0; g < arguments->num_matrices; g++) {
             for (int t = 0; t < lines ; t++){
		 Matrix[g][t][0] = 1 - ((h * t)+(start * h));
		 Matrix[g][t][N] = ((h * t) + (start*h));
		//printf("g: %d, i: %d\n", g, t);	
             }
	}
    }
//printf("end initialize\n");
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
                        {// i NEEDS TO BE CHANGED!!!
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
	//printf("end calculate myrank: %d\n", myrank);
}

static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;
  double** Matrix = arguments->Matrix[results->m];
  MPI_Status status;

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
    int myrank, size, from, to;
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;
    
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
    calculate(myrank, size, &arguments, &results, &options);
 
    for (int i = 0; i < myrank; i++){
	from = from + calculate_lines(i,size, options.interlines);
    }
    from = from + 1;

   for (int i = 0; i <= myrank; i++){
	to = to + calculate_lines(i, size, options.interlines);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    DisplayMatrix(&arguments, &results, &options, myrank, size, from, to);
}
