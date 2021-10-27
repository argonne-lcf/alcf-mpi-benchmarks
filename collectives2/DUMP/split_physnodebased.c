#include <mpi.h>
#include <stdio.h>
#include <pmi.h>
#include <stdlib.h>
#include "xctopo.h"

extern int ppn;
extern int node_size, node_rank;
extern MPI_Comm node_comm, aries_comm;

//#define NUM_NODES 3240 // Number of total compute nodes in Theta


int xctopo_get_mycoords(xctopo_t * topo)
{
  FILE * procfile = fopen("/proc/cray_xt/cname","r");
  if (procfile!=NULL) {

    char a, b, c, d;
    int col, row, cage, slot, anode;

    /* format example: c1-0c1s2n1 c3-0c2s15n3 */
    fscanf(procfile, 
           "%c%d-%d%c%d%c%d%c%d", 
           &a, &col, &row, &b, &cage, &c, &slot, &d, &anode);

#ifdef DEBUG
    fprintf(stderr, "coords = (%d,%d,%d,%d,%d) \n", rank, col, row, cage slot, anode);
#endif
    
    topo->col   = col;
    topo->row   = row;
    topo->cage  = cage;
    topo->slot  = slot;
    topo->anode = anode;

    fclose(procfile);

  } else {

    fprintf(stderr, "xctopo_get_mycoords: fopen has failed! \n");
    return 1;

  }

  return 0;
}


int *coords, *tmp;
int rank;
int _src;
xctopo_t topo;

void freeMem()
{
  free(coords);
  free(tmp);
}



/*
  Assigns ppn
*/




int init()
{
  int i, rc, spawned, pmi_rank, pmi_size;
  PMI_BOOL initialized;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  rc = xctopo_get_mycoords(&topo);

  int col   = topo.col;
  int row   = topo.row;
  int cage  = topo.cage;
  int slot  = topo.slot;
  int anode = topo.anode;


  /* Using PMI to determine ppn */
  rc = PMI_Initialized( &initialized );
  if ( initialized != PMI_TRUE ) 
      rc = PMI_Init( &spawned );
  
  PMI_Get_size( &pmi_size );
  PMI_Get_rank( &pmi_rank );
  PMI_Get_clique_size( &ppn );

  return 0;
  
}


/* 
   Returns two communicators, which cut a partition in two halves at coordinate c
   c = 0 - coordinate X
   c = 1 - coordinate Y
   c = 2 - coordinate Z
    
   Returns:
        0, if split was successful
        1, if given incorrect c
*/

#if 0
int split( int c, MPI_Comm *half0, MPI_Comm *half1, int *ppn )
{
  MPIX_Hardware_t hw;
  MPIX_Hardware(&hw);
  *ppn = hw.ppn;
  
  int size, rank;
  int Nlist0, *list0, Nlist1, *list1;
  MPI_Group World_group, Comm0_group, Comm1_group;

  int * coords = (int*) malloc(sizeof(int) * (hw.torus_dimension+1));
  
  /* Make 2 lists of ranks: list0 to exclude from half0,
     list1 to exclude from list 1 */
  MPI_Comm_size( MPI_COMM_WORLD, (int *)&size );
  list0 = (int*) malloc( sizeof( int ) * size );
  list1 = (int*) malloc( sizeof( int ) * size );

  Nlist0 = 0;
  Nlist1 = 0;

  switch ( c )
  {
      case ( 0 ):
          for ( rank = 0; rank < size; rank++ )
          {
              MPIX_Rank2torus(rank, &coords[0]);
              if (coords[0] < hw.Size[0] / 2 ) 
                  list1[ Nlist1++ ] = rank;
              else
                  list0[ Nlist0++ ] = rank;
          }
          break;
      case ( 1 ):
          for ( rank = 0; rank < size; rank++ )
          {
              MPIX_Rank2torus(rank, &coords[0]);
              if (coords[1] < hw.Size[1] / 2 ) 
                  list1[ Nlist1++ ] = rank;
              else
                  list0[ Nlist0++ ] = rank;

          }
          break;
      case ( 2 ):
          for ( rank = 0; rank < size; rank++ )
          {
              MPIX_Rank2torus(rank, &coords[0]);
              if (coords[2] < hw.Size[2] / 2 ) 
                  list1[ Nlist1++ ] = rank;
              else
                  list0[ Nlist0++ ] = rank;
          }
          break;
      case ( 3 ):
          for ( rank = 0; rank < size; rank++ )
          {
              MPIX_Rank2torus(rank, &coords[0]);
              if (coords[3] < hw.Size[3] / 2 ) 
                  list1[ Nlist1++ ] = rank;
              else
                  list0[ Nlist0++ ] = rank;
          }
          break;
      case ( 4 ):
          for ( rank = 0; rank < size; rank++ )
          {
              MPIX_Rank2torus(rank, &coords[0]);
              if (coords[4] < hw.Size[4] / 2 ) 
                  list1[ Nlist1++ ] = rank;
              else
                  list0[ Nlist0++ ] = rank;
          }
          break;
          
     default:
          free( list0 ); 
          free( list1 ); 
          free( coords );
          return 1;
  }    

  MPI_Comm_group( MPI_COMM_WORLD, &World_group );
  MPI_Group_excl( World_group, Nlist0, list0, &Comm0_group );
  MPI_Group_excl( World_group, Nlist1, list1, &Comm1_group );
  MPI_Comm_create( MPI_COMM_WORLD, Comm0_group, half0 );
  MPI_Comm_create( MPI_COMM_WORLD, Comm1_group, half1 );
  MPI_Group_free( &Comm0_group );
  MPI_Group_free( &Comm1_group );

  free( list0 );
  free( list1 );
  free( coords );
  return 0;
}

#endif

#if 0
int split(int c, MPI_Comm *newcomm, int *ppn, int *color)
{
  
  #if 0
  MPIX_Hardware_t hw;
  MPIX_Hardware(&hw);
  *ppn = hw.ppn;
  #endif
 
  int nid, rank, key, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  rc = PMI_Get_nid(rank, &nid);
  if (rc!=PMI_SUCCESS)
    PMI_Abort(rc,"PMI_Get_nid failed");
  printf("rank %d PMI_Get_nid gives nid %d \n", rank, nid); 

#if 1
  int size, rank, key;
  
  *color = 1;
  //key = rank;
  key = nid;

  // c >= 0 and < NUM_NODES/2 
  if (nid < NUM_NODES / 2)
    *color = 0;

  MPI_Comm_split(MPI_COMM_WORLD, *color, key, newcomm);

#endif
  return 0;
}
#endif
 

// Theta 
// works only when run on the whole machine, to make it work for other cases, communicators should be created using Comm_group
int split2(MPI_Comm *newcomm, int *ppn, int *color)
{

  int rc, nid, rank, key, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  rc = PMI_Get_nid(rank, &nid);
  if (rc!=PMI_SUCCESS)
    PMI_Abort(rc,"PMI_Get_nid failed");
  //printf("rank %d PMI_Get_nid gives nid %d \n", rank, nid); 


  //MPI_Group world_group;
  //MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  

  *color = 1;
  key = rank;

  // c >= 0 and < NUM_NODES/2 
  if (nid < NUM_NODES / 2)
    *color = 0;

 MPI_Comm_split(MPI_COMM_WORLD, *color, key, newcomm);

  return 0;
}
 

// Theta 
// works only when run on the whole machine, to make it work for other cases, communicators should be created using Comm_group
int split4( MPI_Comm *newcomm, int *ppn, int *color)
{

  int rc, nid, rank, key, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  rc = PMI_Get_nid(rank, &nid);
  if (rc!=PMI_SUCCESS)
    PMI_Abort(rc,"PMI_Get_nid failed");
  //printf("rank %d PMI_Get_nid gives nid %d \n", rank, nid); 


  //MPI_Group world_group;
  //MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  

  *color = 3;
  key = rank;

  // c >= 0 and < NUM_NODES/2 
  if (nid < NUM_NODES / 4)
    *color = 0;
  else if (nid >= NUM_NODES / 4 && nid < 2*(NUM_NODES/4))
    *color = 1;
  else if (nid >= 2*(NUM_NODES / 4) && nid < 3*(NUM_NODES/4))
    *color = 2;

 MPI_Comm_split(MPI_COMM_WORLD, *color, key, newcomm);
 return 0;
}
