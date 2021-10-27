/* 
   Argonne Leadership Computing Facility benchmark
   Cray XC version
   Aggregate bandwidth 
   Written by Vitali Morozov <morozov@anl.gov>
   Modified by: Sudheer Chunduri <sudheer@anl.gov>
*/    
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define MAXSKIP 10     // skip first tests to warm-up
#define MAXN 256      // repetition rate for a single pair test
#define MAXTASKS 64
#define XNBORS 3   // max number of neighbors, can be 0, 1, 2, or 3: aries chip is shared by 4 nodes  */

//#define LENGTH 262144
#define LENGTH 524288 //message size chosen such a way that it fits in the L2 cache on tile.
//#define LENGTH 2097152 //message size chosen such a way that it fits in the L2 cache on tile.

int ppn;
int layer_size, layer_rank, node_size, node_rank, report_size, report_rank, aries_size, aries_rank;
MPI_Comm layer_comm, node_comm, report_comm, aries_comm;
MPI_Comm comm = MPI_COMM_WORLD;
//MPI_Request req[2 * XNBORS * MAXTASKS];
MPI_Request req[XNBORS * MAXTASKS];

/* Returns Aggregate bandwidth */
inline double measureAggregate()
{
  int tag = 10, i, j, k, w, N, total;
  char *sb, *rb;
  
  double  d = 0.0;
  double t1;

    posix_memalign( (void **)&sb, 16, sizeof( char ) * LENGTH );    
    posix_memalign( (void **)&rb, 16, sizeof( char ) * LENGTH );    

  //  MPI_Barrier (comm);
  // num_targets = layer_size - 1;
  N = layer_size - 1;
  
  MPI_Barrier( layer_comm );
#if 1
  if ( layer_rank == 0 ) 
  {
      /* source part */
      // heating up! //
      for ( k = 0; k < MAXSKIP; k++ )
      {
        for ( i = 0; i < N; i++ )
        {
          MPI_Isend( sb, LENGTH, MPI_BYTE, i+1, tag, layer_comm, &req[i] );
          //MPI_Irecv( rb, LENGTH, MPI_BYTE, i+1, tag, layer_comm, &req[i+N] );
        }

        //MPI_Waitall( 2*N, req, MPI_STATUSES_IGNORE );
        MPI_Waitall( N, req, MPI_STATUSES_IGNORE );
      }

      t1 = MPI_Wtime();

      for ( k = 0; k < MAXN; k++ )
      {
        for ( i = 0; i < N; i++ )
        {
          MPI_Isend( sb, LENGTH, MPI_BYTE, i+1, tag,    layer_comm, &req[i] );
          //MPI_Irecv( rb, LENGTH, MPI_BYTE, i+1, tag, layer_comm, &req[i+N] );
        }
        //MPI_Waitall( 2*N, req, MPI_STATUSES_IGNORE );
        MPI_Waitall( N, req, MPI_STATUSES_IGNORE );
      }

      t1 = MPI_Wtime() - t1;


      // in MB/s, 1MB = 1e6 B 
      //d = ( (2.0 * (double)LENGTH * (double)N * (double)MAXN ) / (t1) );
      d = ( (1.0 * (double)LENGTH * (double)N * (double)MAXN ) / (t1) );
      d /= 1e9;

    }


  else
  {
      /* receiver part */
     for ( k = 0; k < MAXN+MAXSKIP; k++ )
      {
          total = 0;
          //MPI_Isend( sb, LENGTH, MPI_BYTE, 0, tag, layer_comm, &req[total++] );
          MPI_Irecv( rb, LENGTH, MPI_BYTE, 0, tag,    layer_comm, &req[total++] );

        //MPI_Waitall( 2, req, MPI_STATUSES_IGNORE );
        MPI_Waitall( 1, req, MPI_STATUSES_IGNORE );

      }

  }
#endif
    free( sb );
    free( rb );

  return d;

}



int main( int argc, char *argv[] )
{
  int taskid, ntasks, provided, rc, color;
  double d;
  double maxd, mind, sumd;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  MPI_Comm_size( MPI_COMM_WORLD, &ntasks );

  rc = init();

  d = 0.;
  d = measureAggregate();

  MPI_Reduce( &d, &sumd, 1, MPI_DOUBLE, MPI_SUM, 0, node_comm );

  color = 5;
  if ( layer_rank == 0 && node_rank == 0 ) color = layer_size;
  MPI_Comm_split( MPI_COMM_WORLD, color, taskid, &report_comm );
  MPI_Comm_size( report_comm, &report_size );
  MPI_Comm_rank( report_comm, &report_rank );

  MPI_Reduce( &sumd, &mind, 1, MPI_DOUBLE, MPI_MIN, 0, report_comm );
  MPI_Reduce( &sumd, &maxd, 1, MPI_DOUBLE, MPI_MAX, 0, report_comm );

/*
  if ( layer_rank == 0 )
     printf( "Rank %d: source, time = %lf, aggregate - %lf\n", taskid, t, d);
  if ( layer_rank == 0 && node_rank == 0 )
     printf( "Node Rank %d: source, time = %lf, total aggregate - %lf\n", taskid, t, sumd );
*/

  if ( report_rank == 0 && color < 5 )
     printf( "ppn = %d: Aries %d: nodes %d: min/max aggregate  %lf / %lf GB/s\n", ppn, layer_size, report_size, mind, maxd );
      //printf( "Source rank %d gets %18.12lf GB/s\n", taskid, d );



  MPI_Finalize();
  

  return 0 ;
}
