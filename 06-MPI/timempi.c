#include <mpi.h>
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>
int main (int argc, char **argv)
{
int myrank, size, status;
char buffer[256];
char hostname[16];
struct timeval timeofday;
time_t t;
struct tm * ts;
status = 0;
t = time(NULL);
ts = localtime(&t);

/* MPI Initialisierung, Rangabfrage und Anzahlabfrage */
MPI_Init( &argc, &argv );
MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
MPI_Comm_size( MPI_COMM_WORLD, &size);

/* Ausgabe der einzelnen Prozesse */
if (myrank != 0)
{
	gethostname(hostname, sizeof hostname);
	gettimeofday(&timeofday, NULL);
	snprintf(buffer,sizeof(buffer), "%s: %d-%d-%d %d:%d:%d.%ld\n", hostname,
		 ts->tm_year+1900, ts->tm_mon+1, ts->tm_mday, ts->tm_hour, ts->tm_min,
		 ts->tm_sec, timeofday.tv_usec);
	MPI_Send(&buffer, sizeof(buffer), MPI_CHAR, 0, 99 ,  MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Rang %d beendet jetzt!\n", myrank);
}
else
{
	for (int i = 1; i < size; i++)
	{
	MPI_Recv(&buffer, sizeof(buffer), MPI_CHAR, i, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	printf("%s",buffer);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Rang %d beendet jetzt!\n", myrank);
}
/* Beenden der MPI Routine */
MPI_Finalize();
return 0;
}
