#!/bin/bash
set -euo pipefail

readonly TARGET="filled.gro"
readonly N=100

function runit(){
    readonly GMXNAME=model-test
    export GMX_MAXBACKUP=-1
    grompp -f nve.mdp -c ${TARGET} -p model.top -n qm.ndx -o ${GMXNAME}.tpr -maxwarn 1
    mdrun -v -deffnm ${GMXNAME}
}

function recordit(){
    k=$1
    tail
}

function sample(){
    k=$1
    buildit $k
    for i in $(seq $N)
    do
        runit
        recordit k
    done
}


for k in $(seq 5 30 -w)
do
    sample $k
done
