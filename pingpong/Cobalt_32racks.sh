#!/bin/bash -x
export L1P_POLICY=std
export BG_THREADLAYOUT=1   # 1 - default next core first; 2 - my core first

export PROG=pingpong
export NODES=32768
export RANKS_PER_NODE=2

NPROCS=$((NODES*RANKS_PER_NODE)) 
OUTPUT=N${NODES}_R${RANKS_PER_NODE}_mmps_${NPROCS}
VARS="PAMID_VERBOSE=1:BG_SHAREDMEMSIZE=64:PAMID_EAGER=20000"
rm -f core.* ${OUTPUT}.cobaltlog ${OUTPUT}.output ${OUTPUT}.error
qsub -A Acceptance -n ${NODES} --mode c${RANKS_PER_NODE} -t 0:10:00 -O $OUTPUT --env ${VARS} $PROG

