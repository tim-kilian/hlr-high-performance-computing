#include <unistd.h>
#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>
#include <time.h>


int main (int argc, char **argv)
	{
		int rank, size;
		struct timeval tv;
		struct tm* ptm;
		char time_string[40];
		char hostname[30];
		long mikroseconds, reduce;
		char buffer[100];




		MPI_Init( &argc, &argv );
		MPI_Comm_rank( MPI_COMM_WORLD, &rank );
		MPI_Comm_size( MPI_COMM_WORLD, &size );


	if (rank != 0)
	{
		/* Ausgabe */
		gettimeofday(&tv, NULL);
		gethostname(hostname, sizeof(hostname));
		ptm = localtime (&tv.tv_sec);
		strftime (time_string, sizeof(time_string), "%Y-%m-%d %H:%M:%S",ptm);
		mikroseconds = tv.tv_usec;
		snprintf(buffer, sizeof(buffer), "%s: %s.%06ld\n", hostname, time_string, mikroseconds);
		MPI_Send(&buffer, sizeof(buffer), MPI_CHAR, 0, 999, MPI_COMM_WORLD);
		/*MPI_Send(&mikroseconds, sizeof(mikroseconds), MPI_LONG, 0, 999, MPI_COMM_WORLD);*/
		/*MPI_Barrier(MPI_COMM_WORLD);*/
		/*printf("Rang %d beendet jetzt!\n", rank);*/
	}
	else
	{
		mikroseconds = 999999;
	for(int i = 1; i < size; i++) 
	{
		MPI_Recv(&buffer, sizeof(buffer), MPI_CHAR, i, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		/*MPI_Recv(&mikroseconds, sizeof(mikroseconds), MPI_LONG, i, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);*/
  		printf("%s\n", buffer);
	}	
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&mikroseconds, &reduce, sizeof(mikroseconds), MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
if (rank == 0)
{
	printf("%ld\n", reduce);
}
	/*MPI_Barrier(MPI_COMM_WORLD);*/
	printf("Rang %d beendet jetzt!\n", rank);


MPI_Finalize();
return 0;

	}