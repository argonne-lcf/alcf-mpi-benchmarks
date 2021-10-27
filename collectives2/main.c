/* 
   Argonne Leadership Computing Facility benchmark
   Cray XC version
   Collectives, delivered latency
   Written by Vitali Morozov <morozov@anl.gov>
   Modified by Sudheer Chunduri <sudheer@anl.gov>

   Measures barrier, broadcast, and allreduce delivered latencies
   for MPI_COMM_WORLD, a copy of MPI_COMM_WORLD, and MPI_COMM_WORLD-1
*/    
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

double start, end;


//#define MAXN 500       // repetition rate for a single pair test
#define MAXN 15000       // repetition rate for a single pair test
#define MAXSKIP 20       // skip first tests to warm-up
#define EXCLUDE_RANK 0   // the rank to exclude from World for World-1 benchmark

#define BSIZE 1024*8     // largest broadcast message size, bytes
#define RSIZE 256*8      // largest allreduce buffer size, bytes

void do_barriers ( MPI_Comm comm, double *dt );
void do_broadcast( MPI_Comm comm, int Length, double *buf, double *dt ); 
void do_allreduce( MPI_Comm comm, int Length, double *sbuf, double *rbuf, double *dt ); 

//int split2( MPI_Comm *newcomm, int * color );
//int split4( MPI_Comm *newcomm, int * color );

struct benchmark
{
  const char *comm;
  double time[10];
};
int ppn;

void* _ALLOC_MAIN_ (size_t size)
{
  void* p_buf;
  posix_memalign( (void **)&p_buf, 64, sizeof( char ) * size);
  if (!p_buf)
  {
    printf("ERROR:  Allocating memory %ld bytes failed\n", size);
    fflush(stdout);
    exit(1);
  }
  memset(p_buf, 0, size);
  return p_buf;
}


main( int argc, char *argv[] )
{
  MPI_Comm newcomm, comm;
  int rc, i, j, total_comms, color, provided;
  int key, rank, taskid, ntasks, wsize;
  double *bbuf, *sbuf, *rbuf;
  
  bbuf = _ALLOC_MAIN_(BSIZE + 2*RSIZE);
  sbuf = (double*) ((char*)bbuf + BSIZE);
  rbuf = (double*) ((char*)sbuf + RSIZE);

  static struct benchmark b[] = 
    {
      { "MPI_COMM_WORLD", 0., 0. },
      { "Copy of MPI_COMM_WORLD", 0., 0. },
      { "2Half0 of WORLD", 0., 0. },
      { "2Half1 of WORLD", 0., 0. },
      { "4Half0 of WORLD", 0., 0. },
      { "4Half1 of WORLD", 0., 0. },
      { "4Half2 of WORLD", 0., 0. },
      { "4Half3 of WORLD", 0., 0. },
    };
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &wsize );

  // call to get ppn
  rc = init();

  MPI_Comm_dup(MPI_COMM_WORLD, &newcomm); 
  int comm_size; 

  // comm_world and dup first
  for (i = 0; i < 2; i++)
  {
    if (i == 0)
      comm = MPI_COMM_WORLD;
    else
      comm = newcomm;
    MPI_Comm_size(comm, &comm_size);

    MPI_Barrier(MPI_COMM_WORLD);
    do_barriers (comm, &b[i].time[0]);
    do_broadcast(comm,    2, bbuf, &b[i].time[2]);
    do_broadcast(comm, 1024, bbuf, &b[i].time[4]);
    do_allreduce(comm,    1, sbuf, rbuf, &b[i].time[6]);
    do_allreduce(comm,  256, sbuf, rbuf, &b[i].time[8]);
    
    MPI_Comm_rank(comm, &taskid );
    
    if ( taskid == 0 )
    {
      printf("PPN: %d size: %d %20s: Barrier                 min, max, us: %18.2lf   %18.2lf\n",
             ppn, comm_size, b[i].comm, b[i].time[0], b[i].time[1]);
      printf("PPN: %d size: %d %20s: Broadcast    2 doubles, min, max, us: %18.2lf   %18.2lf\n",
             ppn, comm_size, b[i].comm, b[i].time[2], b[i].time[3]);
      printf("PPN: %d size: %d %20s: Broadcast 1024 doubles, min, max, us: %18.2lf   %18.2lf\n",
             ppn, comm_size, b[i].comm, b[i].time[4], b[i].time[5]);
      printf("PPN: %d size: %d %20s: Allreduce    1 double , min, max, us: %18.2lf   %18.2lf\n",
             ppn, comm_size, b[i].comm, b[i].time[6], b[i].time[7]);
      printf("PPN: %d size: %d %20s: Allreduce  256 double , min, max, us: %18.2lf   %18.2lf\n",
             ppn, comm_size, b[i].comm, b[i].time[8], b[i].time[9]);
      printf("\n");
    }
  }

  MPI_Comm_free(&newcomm);
  total_comms = 2;
  
  // now split into half
    if (rank < wsize / 2)
    color = 0;
    else color = 1;
    key = rank;
    total_comms += color;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &newcomm);
    MPI_Comm_size(newcomm, &comm_size);
    
      MPI_Barrier(MPI_COMM_WORLD);

        do_barriers (newcomm, &b[total_comms].time[0]);
        do_broadcast(newcomm,    2, bbuf, &b[total_comms].time[2]);
        do_broadcast(newcomm, 1024, bbuf, &b[total_comms].time[4]);
        do_allreduce(newcomm,    1, sbuf, rbuf, &b[total_comms].time[6]);
        do_allreduce(newcomm,  256, sbuf, rbuf, &b[total_comms].time[8]);
        
        MPI_Comm_rank(newcomm, &taskid);
        
        if ( taskid == 0 )
        {
          printf("PPN: %d size: %d %20s: Barrier                 min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[0],
                 b[total_comms].time[1]);
          printf("PPN: %d size: %d %20s: Broadcast    2 doubles, min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[2],
                 b[total_comms].time[3]);
          printf("PPN: %d  size: %d %20s: Broadcast 1024 doubles, min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[4],
                 b[total_comms].time[5]);
          printf("PPN: %d  size: %d %20s: Allreduce    1 double , min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[6],
                 b[total_comms].time[7]);
          printf("PPN: %d  size: %d %20s: Allreduce  256 double , min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[8],
                 b[total_comms].time[9]);
          printf("\n");
        }
    MPI_Comm_free(&newcomm);


  total_comms = 4;
  // now split into quater
     key = rank;

  if (rank < wsize / 4)
    color = 0;
  else if ((rank >= (wsize / 4)) && (rank < 2*(wsize/4)))
    color = 1;
  else if ((rank >= 2*(wsize / 4)) && (rank < 3*(wsize/4)))
    color = 2;
  else if ((rank >= 3*(wsize/4)) && (rank < wsize) )
    color = 3;

    total_comms += color;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &newcomm);
    MPI_Comm_size(newcomm, &comm_size);

      MPI_Barrier(MPI_COMM_WORLD);




        do_barriers (newcomm, &b[total_comms].time[0]);
        do_broadcast(newcomm,    2, bbuf, &b[total_comms].time[2]);
        do_broadcast(newcomm, 1024, bbuf, &b[total_comms].time[4]);
        do_allreduce(newcomm,    1, sbuf, rbuf, &b[total_comms].time[6]);
        do_allreduce(newcomm,  256, sbuf, rbuf, &b[total_comms].time[8]);
        
        MPI_Comm_rank(newcomm, &taskid);
        
        if ( taskid == 0 )
        {
          printf("PPN: %d  size: %d %20s: Barrier                 min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[0],
                 b[total_comms].time[1]);
          printf("PPN: %d  size: %d %20s: Broadcast    2 doubles, min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[2],
                 b[total_comms].time[3]);
          printf("PPN: %d  size: %d %20s: Broadcast 1024 doubles, min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[4],
                 b[total_comms].time[5]);
          printf("PPN: %d  size: %d %20s: Allreduce    1 double , min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[6],
                 b[total_comms].time[7]);
          printf("PPN: %d  size: %d %20s: Allreduce  256 double , min, max, us: %18.2lf   %18.2lf\n",
                 ppn, comm_size, b[total_comms].comm, b[total_comms].time[8],
                 b[total_comms].time[9]);
          printf("\n");
        }

    MPI_Comm_free(&newcomm);

  free(bbuf);
  MPI_Finalize();
  return 0;
}

void do_barriers( MPI_Comm comm, double *dt )    
{
  int k;
  double d;
  start = end = 0;
  
  for ( k = 0; k < MAXSKIP; k++ ) MPI_Barrier( comm );
    
  start = MPI_Wtime();
  for ( k = 0; k < MAXN; k++ )
  {
    MPI_Barrier( comm );
  }
  end = MPI_Wtime() - start;
  d = 1e6* end / (double)MAXN;
    
  MPI_Reduce( &d, &dt[0], 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &d, &dt[1], 1, MPI_DOUBLE, MPI_MAX, 0, comm );

  return;
}

void do_broadcast( MPI_Comm comm, int Length, double *buf, double *dt )    
{
  int k;
  double d;
  start = end = 0;
  MPI_Barrier( comm );
    
  for ( k = 0; k < MAXSKIP; k++ )
    MPI_Bcast( buf, Length, MPI_DOUBLE_PRECISION, 0, comm ); 
    
  start = MPI_Wtime();
  for ( k = 0; k < MAXN; k++ )
  {
    MPI_Bcast( buf, Length, MPI_DOUBLE, 0, comm );
  }
  end = MPI_Wtime() - start;

  d = 1e6 * end / (double)MAXN;
    
  MPI_Reduce( &d, &dt[0], 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &d, &dt[1], 1, MPI_DOUBLE, MPI_MAX, 0, comm );

  return;
}

void do_allreduce( MPI_Comm comm, int Length, double *sbuf, double *rbuf, double *dt )    
{
  int k;
  double d;
  start = end = 0;
  MPI_Barrier( comm );
    
  for ( k = 0; k < MAXSKIP; k++ )
    MPI_Allreduce( sbuf, rbuf, Length, MPI_DOUBLE, MPI_SUM, comm ); 
    
  start = MPI_Wtime();
  for ( k = 0; k < MAXN; k++ )
  {
    MPI_Allreduce( sbuf, rbuf, Length, MPI_DOUBLE, MPI_SUM, comm );
  }
  end = MPI_Wtime() - start;
  d = 1e6* end / (double)MAXN;
    
  MPI_Reduce( &d, &dt[0], 1, MPI_DOUBLE, MPI_MIN, 0, comm );
  MPI_Reduce( &d, &dt[1], 1, MPI_DOUBLE, MPI_MAX, 0, comm );

  return;
}
