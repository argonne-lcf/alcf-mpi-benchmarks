/* 
   Argonne Leadership Computing Facility benchmark
   Cray XC version
   Aggregate bandwidth 
   Written by Vitali Morozov <morozov@anl.gov>
*/    
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define MAXN 100       // repetition rate for a single pair test
#define MAXSKIP 10     // skip first tests to warm-up


//#define N 3         // number of neighbors
//#define RPN 4        // number of ranks per node

#define LENGTH 524288
int ppn;
int layer_size, layer_rank, node_size, node_rank, report_size, report_rank, aries_size, aries_rank;
MPI_Comm layer_comm, node_comm, report_comm, aries_comm;
MPI_Comm comm = MPI_COMM_WORLD;

int getranks(int *, int *);

int main( int argc, char *argv[] )
{
  int rc, provided, N, taskid, ntasks, is1, i, L, k, Receiver, pmi_size, nodes;
  int *ir1;
  char *sb, *rb;
  double d, d1;
  MPI_Status *stat;
  MPI_Request *req;
  double t1;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  MPI_Comm_size( MPI_COMM_WORLD, &ntasks );

  rc = init();

  #if 0
  PMI_Get_size( &pmi_size );
  PMI_Get_clique_size( &ppn );

  nodes = ntasks / ppn;
  N = nodes - 1;
  if ( N < 1 ) 
  {
      printf( "Increase the number of nodes: running on %d nodes now\n", nodes );
      MPI_Finalize();
      return 0;
  }    

  ir1 = (int *)malloc( N * sizeof( int ) );
  stat = (MPI_Status *)malloc( 2 * N * sizeof( MPI_Status ) );
  req = (MPI_Request *)malloc( 2 * N * sizeof( MPI_Request ) );

  if ( taskid == 0 ) printf( "MPI size = %d, Nodes = %d, RPN = %d\n", ntasks, nodes, ppn );
  #endif

  /*
  getranks( N, &is1, ir1, &ppn );
  MPI_Finalize();
  return 0;
  */
  
  
  if ( getranks( &is1, ir1 ) == 0 )
  {
    Receiver = 0;
    for ( i = 0; i < N; i++ )
      if ( taskid == ir1[i] )
        Receiver++;

    d = 0.0;

    posix_memalign( (void **)&sb, 16, sizeof( char ) * LENGTH );    
    posix_memalign( (void **)&rb, 16, sizeof( char ) * LENGTH );    

    MPI_Barrier (comm);

    if ( taskid == is1 )
    {
      #if 1
      for ( i = 0; i < N; i++ )
        printf("%d => %d\n", is1, ir1[i] );
      #endif
      // heating up! //
      for ( k = 0; k < MAXSKIP; k++ )
      {   
        for ( i = 0; i < N; i++ )
        {
          MPI_Isend( sb, LENGTH, MPI_BYTE, ir1[i], is1, comm, &req[i] );
          MPI_Irecv( rb, LENGTH, MPI_BYTE, ir1[i], ir1[i], comm, &req[i+N] );
        }
                
        MPI_Waitall( 2*N, req, stat );
      }
            
      t1 = MPI_Wtime();
            
      for ( k = 0; k < MAXN; k++ )
      {
        for ( i = 0; i < N; i++ )
        {
          MPI_Isend( sb, LENGTH, MPI_BYTE, ir1[i], is1,    comm, &req[i] );
          MPI_Irecv( rb, LENGTH, MPI_BYTE, ir1[i], ir1[i], comm, &req[i+N] );
        }
        MPI_Waitall( 2*N, req, stat );
      }
            
      t1 = MPI_Wtime() - t1;


      // in MB/s, 1MB = 1e6 B 
      d = ( (2.0 * (double)LENGTH * (double)N * (double)MAXN ) / (t1) );
      d /= 1e9;

      printf( "Source rank %d gets %18.12lf GB/s\n", taskid, d ); 
    }

    if ( Receiver )
    {
      for ( k = 0; k < MAXN+MAXSKIP; k++ )
      {
        int total = 0;
        for ( i = 0; i < Receiver; i++ )
        {
          MPI_Isend( sb, LENGTH, MPI_BYTE, is1, taskid, comm, &req[total++] );
          MPI_Irecv( rb, LENGTH, MPI_BYTE, is1, is1,    comm, &req[total++] );
        }

        MPI_Waitall( total, req, stat );

      }
    }

    MPI_Reduce( &d, &d1, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
    if ( taskid == 0 )
      printf("PPN: %d Aggregate BW (GB/s): %18.2lf\n", ppn, d1 );

    free( sb );
    free( rb );
  }

  free( ir1 );
  free( stat );
  free( req  );
  MPI_Finalize();
  

  return 0 ;
}
