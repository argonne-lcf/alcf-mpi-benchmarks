#!/bin/bash -x
export L1P_POLICY=std
export BG_THREADLAYOUT=1   # 1 - default next core first; 2 - my core first

export PROG=aggregate
export NODES=512

rm -f core.* 

#for RANKS_PER_NODE in 1 4 8 16 64 
for RANKS_PER_NODE in 2
do

    NPROCS=$((NODES*RANKS_PER_NODE)) 
    OUTPUT=N${NODES}_R${RANKS_PER_NODE}_mmps_${NPROCS}
    VARS="PAMID_VERBOSE=1:BG_SHAREDMEMSIZE=64"
    rm -f ${OUTPUT}.cobaltlog ${OUTPUT}.output ${OUTPUT}.error
    qsub -A Acceptance -n ${NODES} --mode c${RANKS_PER_NODE} -t 0:10:00 -O $OUTPUT --env ${VARS} $PROG
done

