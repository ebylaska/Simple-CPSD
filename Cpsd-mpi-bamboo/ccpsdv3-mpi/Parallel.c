/* Parallel.c
 *
 *  Created on: Dec 9, 2013
 *      Author: bylaska
 */
#include	<stdio.h>
#include	"mpi.h"

static int taskid,np;


int Parallel_np() {return np;}
int Parallel_taskid() {return taskid;}

void Parallel_init(int *argc, char ***argv)
{
   if (MPI_Init(argc,argv)!=MPI_SUCCESS)
      printf("Parallel_init error: MPI_Init failed\n");

   if (MPI_Comm_rank(MPI_COMM_WORLD,&taskid)!=MPI_SUCCESS)
	  printf("Parallel_init error: MPI_Comm_rank failed\n");

   if (MPI_Comm_size(MPI_COMM_WORLD,&np)!=MPI_SUCCESS)
	  printf("Parallel_init error: MPI_Comm_rank failed\n");
}

void Parallel_end()
{
   if (MPI_Finalize()!=MPI_SUCCESS)
      printf("Parallel_end error: MPI_Finalize failed\n");
}

