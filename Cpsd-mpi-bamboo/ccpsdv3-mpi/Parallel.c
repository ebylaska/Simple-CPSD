/* Parallel.c
 *
 *  Created on: Dec 9, 2013
 *      Author: bylaska
 */
#include	<stdio.h>
#include	"mpi.h"
#include	"float.h"

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

void Parallel_ibcast(int psend, int sz, int *ivalue)
{
   if (MPI_Bcast(ivalue,sz,MPI_INTEGER,psend,MPI_COMM_WORLD)!=MPI_SUCCESS)
      printf("Parallel_ibcast error: MPI_Bcast failed\n");
}

void Parallel_rbcast(int psend, int sz, REAL *value)
{
   if (MPI_Bcast(value,sz,MPI_REAL_PRECISION,psend,MPI_COMM_WORLD)!=MPI_SUCCESS)
      printf("Parallel_rbcast error: MPI_Bcast failed\n");
}


REAL Parallel_SumAll(REAL value_in)
{
   REAL value_out=0.0;
   if (MPI_Allreduce(&value_in,&value_out,1,MPI_REAL_PRECISION,MPI_SUM,MPI_COMM_WORLD)!=MPI_SUCCESS)
      printf("Parallel_SumAll error: MPI_Allreduce failed\n");

   return value_out;
}
