/* 
   Argonne Leadership Computing Facility benchmark
   Cray XC version
   Intranode ping-pong
   Written by Vitali Morozov <morozov@anl.gov>
   Modified by: Daniel Faraj <faraja@us.ibm.com>
   Modified by: Sudheer Chunduri <sudheer@anl.gov>

   Version 0.81: Aug 23, 2012. Only sender must measure time!

   (BG specific: Must run on at least 512 nodes with 2 ranks per node)
   (Cray Aries specific: switch off the adaptive routing, use minimum routing to perform controlled experiments.)
*/    
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXL 3            // number of message sizes 
#define MAXN 150000       // Do not change the repetion rate!

int ppn;
int layer_size, layer_rank, node_size, node_rank, report_size, report_rank, aries_size, aries_rank;
int chasiss_size, chasiss_rank, cabinet_size, cabinet_rank, intergroup_size, intergroup_rank, intragroup_size, intragroup_rank;
MPI_Comm chasiss_comm, cabinet_comm, intragroup_comm, intergroup_comm, layer_comm, node_comm, report_comm, aries_comm;
MPI_Comm comm = MPI_COMM_WORLD;


//char tests[3][16] = {"Intranode", "Nearest", "Farthest"};
char *sb, *rb;
double t = 1e9;

double communicate(int taskid, int is1, int ir1, int L)
{
  int k;

  t = 1e9;
  
  MPI_Barrier (MPI_COMM_WORLD);
  
  if (taskid == 0)
  {
    t = MPI_Wtime();
    for (k = 0; k < MAXN; k++)
    {   

      MPI_Send( sb, L, MPI_CHAR, ir1, is1, MPI_COMM_WORLD );
      MPI_Recv( rb, L, MPI_CHAR, ir1, ir1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    t = MPI_Wtime() - t;
    t = t * 1e6 / MAXN / 2.0;
  }

  else if (taskid == 1)
  {
    for (k = 0; k < MAXN; k++)
    { 
      MPI_Recv( rb, L, MPI_CHAR, is1, is1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send( sb, L, MPI_CHAR, is1, ir1, MPI_COMM_WORLD );
    }
  }

  return t;
}

main( int argc, char *argv[] )
{
  int taskid, ntasks, LENGTH[MAXL], is1 = -1, ir1 = -1, i, j, L;
  double t, total;
  
  LENGTH[0] = 0;
  LENGTH[1] = 4096;
  LENGTH[2] = 65536;
  
  posix_memalign( (void **)&sb, 64, sizeof( char ) * LENGTH[ MAXL-1 ] );
  posix_memalign( (void **)&rb, 64, sizeof( char ) * LENGTH[ MAXL-1 ] );

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  MPI_Comm_size( MPI_COMM_WORLD, &ntasks );
  
  if ( taskid == 0 )
    printf("%12s %6s  %10s\n", "Comm type", "msize", "Latency(us)");

  MPI_Barrier(MPI_COMM_WORLD);
   
  //for (j = 0; j < 3; j++)
  {
    //getranks(j, &is1, &ir1 );
    for ( i = 0; i < MAXL; i++ )
    {
      L = LENGTH[ i ];

      // warm up
      communicate(taskid, 0, 1, L);
      
      t = communicate(taskid, 0, 1, L);

  /*    if ( j == 0 )
      {
         // intranode, all nodes are participating and minimal time is reported.
         MPI_Reduce(&t, &total, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

         if ( taskid == 0 )
             printf("%12s %6d: %10.4lf\n", tests[j], L, total);
      }
       else
   */
         if ( taskid == 0 )
             //printf("%12s %6d: %10.4lf\n", tests[j], L, t );
             printf(" %6d: %10.4lf\n", L, t );
    }
  }

  free(sb);
  free(rb);
  
  MPI_Finalize();
  return 0;
}
