/* 
   Argonne Leadership Computing Facility benchmark
   Cray XC version
   Messaging rate 
   Written by Vitali Morozov <morozov@anl.gov>
   Modified by: Daniel Faraj <faraja@us.ibm.com>
*/    
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


#define MAXN 256      // repetition rate for a single pair test
#define MAXTASKS 64
#define MAXWIN 8
#define XNBORS 3   // max number of neighbors, can be 0, 1, 2, or 3: aries chip is shared by 4 nodes  */

MPI_Comm comm = MPI_COMM_WORLD;
MPI_Request req[2 * XNBORS * MAXTASKS * MAXWIN];
int targets[XNBORS * MAXTASKS];
double MMPS_W[ MAXWIN ];

int ppn;
int layer_size, layer_rank, node_size, node_rank, report_size, report_rank, aries_size, aries_rank;
MPI_Comm layer_comm, node_comm, report_comm, aries_comm;

double maxt;
double t;

//inline void measureMMPS( int taskid, int source, int ppn, int num_targets, int receiver, int window)
/* Returns MMPS */
inline double measureMMPS( int window )
{
  int tag = 10, i, j, k, w, num_targets, iter, total;
  char sb, rb;
  double mmps;

  num_targets = layer_size - 1;
  
  MPI_Barrier( layer_comm );

  if ( layer_rank == 0 ) 
  {
      /* source part */
     
      t = MPI_Wtime();
      for ( iter = 0; iter < MAXN; iter ++ )
      {
          total = 0;
	  for ( w = 0; w < window; w++ )
              for ( j = 0; j < num_targets; j++ )
                   MPI_Irecv( &rb, 0, MPI_CHAR, j+1, tag, layer_comm, &req[total++]);
    
          for ( w = 0; w < window; w++ )
	      for ( j = 0; j < num_targets; j++ )
    	           MPI_Isend( &sb, 0, MPI_CHAR, j+1, tag, layer_comm, &req[total++]);
         
          MPI_Waitall( 2 * num_targets * w, req, MPI_STATUSES_IGNORE );
      }
      
      t = MPI_Wtime() - t;

      mmps =  (double)( 2 * num_targets * w * MAXN  ) / t / 1e6;
  }
  else
  {
      /* receiver part */
      for ( iter = 0; iter < MAXN; iter ++ )
      {
          total = 0;

	  for ( w = 0; w < window; w++ )
              MPI_Irecv( &rb, 0, MPI_CHAR, 0, tag, layer_comm, &req[total++]);
          
	  for ( w = 0; w < window; w++ )
	      MPI_Isend( &sb, 0, MPI_CHAR, 0, tag, layer_comm, &req[total++]);
    
          MPI_Waitall( 2 * w, req, MPI_STATUSES_IGNORE );
      }
      mmps = 0.;
  }

  return mmps;

}

int main( int argc, char *argv[] )
{

 // int ppn, num_targets, taskid, ntasks, source, i, j, window, receiver, best_win, provided;
  int taskid, ntasks, provided, source, rc, color, w;
  
  double best_mrate, mrate[MAXWIN+1], max_rate = 0.0, d, maxd, mind, sumd;
  
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank( MPI_COMM_WORLD, &taskid );
  MPI_Comm_size( MPI_COMM_WORLD, &ntasks );
  
  rc = init();

  for ( w = 1; w < MAXWIN; w++ )
      MMPS_W[ w ] = measureMMPS( w );

  d = 0.;
  for ( w = 1; w < MAXWIN; w++ ) 
      if ( d < MMPS_W[ w ] ) d = MMPS_W[ w ];

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
     printf( "Rank %d: source, time = %lf, mmps - %lf\n", taskid, t, d);


  if ( layer_rank == 0 && node_rank == 0 )
     printf( "Node Rank %d: source, time = %lf, total mmps - %lf\n", taskid, t, sumd );
*/

  if ( report_rank == 0 && color < 5 )
     printf( "ppn = %d: Aries %d: nodes %d: min/max mmps  %lf / %lf\n", ppn, layer_size, report_size, mind, maxd );
     




  MPI_Finalize();
  return 0;

/*
  if (taskid == 0)
  {
    printf("PPN   XNBORS   WIN   MMPS\n");
    printf("-------------------------\n");
  }
  
  




  //if (ppn == 64) iter = 64;


  for (i = 9; i <= XNBORS; i++)
  {
    num_targets = i * ppn;
    getranks(i, targets);
    
    for (receiver = 0, j = 0; j < num_targets; j++)
      if (taskid == targets[j]) receiver++;

    for (window = 1; window <= MAXWIN; window++)
    {
      // warmup
      measureMMPS(taskid, source, ppn, num_targets, receiver, window);

      measureMMPS(taskid, source, ppn, num_targets, receiver, window);

      // measurement
      mrate[window] = (double)
        (2 * ppn * iter * window * num_targets) / maxt / 1e6;
    }

    best_mrate = mrate[1];
    best_win = 1;
    for (j = 2; j < MAXWIN; j++)
      if (best_mrate < mrate[j])
      {
        best_mrate = mrate[j];
        best_win = j;
      }

    if (taskid == 0)
      printf("%3d   %6d   %3d   %.2f\n", ppn, i, best_win, best_mrate);
    
    if (best_mrate > max_rate)
    {
      max_rate = best_mrate;
      best_mrate *= 1.0;
    }    
  }
  if (taskid == 0)
    printf("Max Mrate for ppn %3d: %.2f\n", ppn, max_rate);
  

  freeMem();

  MPI_Finalize();
*/

  return 0 ;
}
