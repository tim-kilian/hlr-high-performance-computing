#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int*
init (int N)
{
  // Todo
  int* buf = malloc(sizeof(int) * N);

  srand(time(NULL));

  for (int i = 0; i < N; i++)
  {
    // Do not modify % 25
    buf[i] = rand() % 25;
  }

  return buf;
}

int*
circle (int* buf)
{
  // Todo
  return buf;
}

int
main (int argc, char** argv)
{
  char arg[256];
  int N;
  int rank;
  int* buf;

  if (argc < 2)
  {
    printf("Arguments error!\n");
    return EXIT_FAILURE;
  }

  sscanf(argv[1], "%s", arg);

  // Array length
  N = atoi(arg);
  buf = init(N);

  // Todo: myrank
  rank = 0;

  printf("\nBEFORE\n");

  for (int i = 0; i < N; i++)
  {
    printf ("rank %d: %d\n", rank, buf[i]);
  }

  circle(buf);

  printf("\nAFTER\n");

  for (int j = 0; j < N; j++)
  {
    printf ("rank %d: %d\n", rank, buf[j]);
  }

  return EXIT_SUCCESS;
}
