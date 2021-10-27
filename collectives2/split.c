#include <mpi.h>
#include <stdio.h>
#include <pmi.h>
#include <stdlib.h>
#include "xctopo.h"

extern int ppn;
extern int node_size, node_rank;
extern MPI_Comm node_comm, aries_comm;


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



 
