#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

#define RECEIVER 0
#define LEN 256

char buffer[LEN];
int size, rank;

int status;
int microseconds = 999999;

void receiveMessage() {
    for (int i = 1; i < size; i++) {
        MPI_Recv(buffer, LEN, MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("%s\n", buffer);
    }
}

const char* hostname() {
    static char hostname[32];
    gethostname(hostname, 32);
    return hostname;
}

const char* timestamp() {
    struct timeval timeval;
    gettimeofday(&timeval, NULL);

    time_t rawtime = timeval.tv_sec;

    struct tm* timeinfo = localtime(&rawtime);

    char timebuffer[32];
    strftime(timebuffer, 32, "%F %T", timeinfo);

    microseconds = timeval.tv_usec;

    static char timestamp[32];
    snprintf(timestamp, 32, "%s.%06ld", timebuffer, timeval.tv_usec);
    return timestamp;
}

void sendMessage() {
    snprintf(buffer, LEN, "%s: %s", hostname(), timestamp());
    // Send the message to the receiver process
    MPI_Send(buffer, LEN, MPI_BYTE, RECEIVER, 0, MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
    // MPI Initialisierung
    MPI_Init(&argc, &argv);
    // MPI Anzahlabfrage: size: Anzahl der Prozesse
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // MPI Rangabfrage: Rand des Prozesses
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Ausgabe der einzelnen Prozesse
    if (rank == RECEIVER) {
        receiveMessage();
    } else {
        sendMessage();
    }

    MPI_Reduce(&microseconds, &status, 1, MPI_INT, MPI_MIN, RECEIVER, MPI_COMM_WORLD);

    if (rank == RECEIVER) {
        printf("%d\n", status);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Rang %d beendet jetzt!\n", rank);

    // Beenden der MPI Routine
    MPI_Finalize();
    return 0;
}
