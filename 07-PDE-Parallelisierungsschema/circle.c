#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int*
init (int N)
{
  // Todo
  int* buf = malloc(sizeof(int) * N);

  srand(time(NULL));
    // Anzahl der Prozesse viele Randomvariablen werden auf buf geschrieben
  for (int i = 0; i < N; i++)
  {
    // Do not modify % 25
    buf[i] = rand() % 25;
  }

  return buf;
}


int*
displ (int N,int size)
{
  int* displacement = malloc(sizeof(int) * size);
   
  for (int i = 0; i < size; i++)
  {
      displacement[i] = (i) * (N/size) ;
  }
   return displacement;
}

int*
sendc (int N,int size)
{
   int* sendcounts = malloc(sizeof(int) * size);
   int rest;
   rest = N % size;
   for (int i = 0; i < size; i++)
   { 
      if (rest == 0)
      {
          sendcounts[i] = (N/size);
      }
      else 
      {
          if (i == (N-1))
          {
	      sendcounts[i] = rest;
          }
          else
          {
	      sendcounts[i] = (N/size);
          }
      } 
   }
   return sendcounts;
}


// Routine gibt Eintraege aus buf weiter
int*
circle (int myrank, int size , int* chunk)
{
  //int* recv = malloc(sizeof(chunk));
  int* first_value = malloc(sizeof(int));
  int abbruch, sum;
  abbruch = 0;
  if (myrank == 0)
  {
    *first_value = chunk[0];
    MPI_Send(first_value, sizeof(int), MPI_INT, size-1, 42, MPI_COMM_WORLD);
  }
  if (myrank == size-1)
  {
    MPI_Recv(first_value, sizeof(int), MPI_INT, 0, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }    
  // Kreisroutine bauen
  for (int i = 0; i < size; i++)
  {
    if (myrank == size-1)
    {
       MPI_Send(chunk, sizeof(chunk), MPI_INT, 0 , 43, MPI_COMM_WORLD);
    }
    else
    {
      MPI_Send(chunk, sizeof(chunk), MPI_INT, myrank + 1, 43, MPI_COMM_WORLD);
    }
    if (myrank == 0)
    {
       MPI_Recv(chunk, sizeof(chunk), MPI_INT, size -1, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
       MPI_Recv(chunk, sizeof(chunk), MPI_INT, myrank-1, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == size-1)
    {  
       if (*first_value == chunk[0])
       {
          abbruch = 1;
       } 
    }
    MPI_Allreduce(&abbruch, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sum == 1) {printf("SUMME %d\n", sum); break;}
  }
  return chunk;
}

int
main (int argc, char** argv)
{
  char arg[256];
  int N;
  int rank, size;
  int* buf;
  
  if (argc < 2)
  {
    printf("Arguments error!\n");
    return EXIT_FAILURE;
  }
// Einlesen der Anzahl an Prozessen aus Kommandozeile
  sscanf(argv[1], "%s", arg);

  // Array length
  N = atoi(arg);
  buf = init(N);
  // Initialisieren des MPI Prozesses
  MPI_Init( &argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  for (int j= 0; j < size; j++) {
  	if (j == rank) {
  		printf("\nBEFORE\n");

  		for (int i = 0; i < N; i++)
  		{
    			printf ("rank %d: %d\n", rank, buf[i]);
  		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
  }

  int* sendcounts;
  int* displacement;
  
  displacement = displ(N, size);
  sendcounts = sendc(N, size);

  int* chunk=malloc(sizeof(int) * sendcounts[rank]);
  MPI_Scatterv(buf, sendcounts, displacement, MPI_INT, chunk, sendcounts[rank+1], MPI_INT,0, MPI_COMM_WORLD);

  circle(rank,size, chunk);
  //MPI_Barrier(MPI_COMM_WORLD);

  //MPI_Gatherv(chunk, sendcounts[rank+1], MPI_INT, buf, sendcounts, displacement, MPI_INT , 0, MPI_COMM_WORLD);
  for (int i = 0; i < size; i++)
  {  
     MPI_Gatherv(chunk, sendcounts[rank], MPI_INT, buf, sendcounts, displacement, MPI_INT, i, MPI_COMM_WORLD);
  }
  for (int j= 0; j < size; j++) {
        if (j == rank) {
                printf("\nAFTER\n");

                for (int i = 0; i < N; i++)
                {
                        printf ("rank %d: %d\n", rank, buf[i]);
                }
        }
	MPI_Barrier(MPI_COMM_WORLD);
  }

  //MPI_Finalize();
  return EXIT_SUCCESS;
}
