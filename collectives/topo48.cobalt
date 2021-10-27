#!/bin/bash -x
#COBALT -A Performance
#COBALT -n 48
#COBALT -t 1:00
#COBALT -O topo48huge
#COBALT -q cache-quad
##COBALT --attrs location=1236-1284

module load craype-hugepages8M

rpn=1
aprun -N $rpn -n $((COBALT_JOBSIZE*rpn)) ./collectives

#rpn=2
#aprun -N $rpn -n $((COBALT_JOBSIZE*rpn)) ./collectives

rpn=4
aprun -N $rpn -n $((COBALT_JOBSIZE*rpn)) ./collectives

#rpn=8
#aprun -N $rpn -n $((COBALT_JOBSIZE*rpn)) ./collectives

rpn=16
aprun -N $rpn -n $((COBALT_JOBSIZE*rpn)) ./collectives

#rpn=32
#aprun -N $rpn -n $((COBALT_JOBSIZE*rpn)) ./collectives

rpn=64
aprun -N $rpn -n $((COBALT_JOBSIZE*rpn)) ./collectives



exit $?

