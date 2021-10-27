#include <mpi.h>
#include <stdio.h>
#include <pmi.h>
#include <stdlib.h>
#include "xctopo.h"

extern int ppn;
extern int layer_size, layer_rank, node_size, node_rank, aries_size, aries_rank;
extern MPI_Comm layer_comm, node_comm, aries_comm;




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
  Returns the source task and the list of the nearest target tasks (1 hop distance) 
  The source task should have 4 nodes connected to the same Aries chip and located on node 0  
  The target task should have 4 nodes connected to the same Aries chip and located on node 1, 2, or 3

  The layer_comm communicator entirely described this scenario
  

  Returns:
  0, if source and target are assigned
  1, if source node does not contain tasks
  2, if one of the target nodes does not contain tasks
*/




int init()
{
  int i, rc, spawned, pmi_rank, pmi_size, node_enc;
  char bin_buf[33];
  PMI_BOOL initialized;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  rc = xctopo_get_mycoords(&topo);

  int col   = topo.col;
  int row   = topo.row;
  int cage  = topo.cage;
  int slot  = topo.slot;
  int anode = topo.anode;

  coords = (int*) malloc( sizeof(int) * 6 ); // column, row, cage, slot, anode, one more for cpu on a node
  tmp = (int*) malloc(sizeof(int) * 3);     // 4 nodes are located on an Aries, 3 nearest neighbors

  /* Using PMI to determine ppn */
  rc = PMI_Initialized( &initialized );
  if ( initialized != PMI_TRUE ) 
      rc = PMI_Init( &spawned );
  
  PMI_Get_size( &pmi_size );
  PMI_Get_rank( &pmi_rank );
  PMI_Get_clique_size( &ppn );


  /* 
  Create a single integer bit encoded value:
  column: 0-11, 4 bits
  rows:   0- 1, 1 bit
  cage:   0- 3, 2 bits
  slot:   0-15, 4 bits
  node:   0-3 , 2 bits
  
  ccccrggssssnn
  */
  node_enc = 0;
  node_enc = node_enc | col;
  node_enc = (node_enc << 1) | row;
  node_enc = (node_enc << 2) | cage;
  node_enc = (node_enc << 4) | slot;
  node_enc = (node_enc << 2) | anode; 

  /* node_comm include all processes on the same node */
  MPI_Comm_split( MPI_COMM_WORLD, node_enc, rank, &node_comm );
  MPI_Comm_size( node_comm, &node_size );
  MPI_Comm_rank( node_comm, &node_rank );

  /* aries_com includes all processes connected to the same aries chip */
  //MPI_Comm_split( MPI_COMM_WORLD, node_enc | 0x3, rank, &aries_comm );
  MPI_Comm_split( MPI_COMM_WORLD, node_enc & 0x1FFC, rank, &aries_comm );
  MPI_Comm_size( aries_comm, &aries_size );
  MPI_Comm_rank( aries_comm, &aries_rank );

  /* layer_comm includes all processes from aries communicator with the same node_rank */
  MPI_Comm_split( aries_comm, node_rank, rank, &layer_comm );
  MPI_Comm_size( layer_comm, &layer_size );
  MPI_Comm_rank( layer_comm, &layer_rank );

  //printf("%d: topology coords = (%d,%d,%d,%d,%d): ppn = %d, aries %d / %d, node %d / %d\n", rank, col, row, cage, slot, anode, ppn, aries_rank, aries_size, node_rank, node_size );
 
/* 
  if ( layer_rank == 0 )
      printf("%d: topology coords = (%d,%d,%d,%d,%d): ppn = %d, source within %d processes\n", rank, col, row, cage, slot, anode, ppn, layer_size );
  else
      printf("%d: topology coords = (%d,%d,%d,%d,%d): ppn = %d, target within %d processes\n", rank, col, row, cage, slot, anode, ppn, layer_size );
*/

  /* source is always rank 0 of communicator layer_comm, no need to return */
  
  if ( layer_size == 1 ) return 1; /* the aries communicator does not have enough nodes to run the benchmark */

  return 0;
  
}

#if 0
/* 
  With benchmark running on layer_comm, the number of neighbors is "layer_size - 1", and targets[] are 1 ... layer_size - 1
*/
void getranks( int *XNbors, int *targets)
{
  int i, j, k = 0;
  int rank;

  *XNbors = layer_size - 1;

  for (i = 0; i < *XNbors; i++)
      targets[i] = i + 1;
  
  return;
}
#endif
