#!/bin/bash
set -euo pipefail

readonly base=$(readlink -f $(dirname $0))

readonly TEMPLATE="${base}/template.gro"


readonly TARGET="filled.gro"

readonly output="${base}/run.dat"
readonly states="cs.dat"


readonly N=4000
readonly int=0.2

readonly GMXNAME=model-test
function runit(){
    rm ${states} || true
    export GMX_MAXBACKUP=-1
    grompp -f nve.mdp -c ${TARGET} -p model.top -n qm.ndx -o ${GMXNAME}.tpr -maxwarn 1
    mdrun -v -deffnm ${GMXNAME}
}

function buildit(){
    k=$1
    vel=$(printf "%8.4f" $(echo "$k * 1.093845650" | bc -l))
    sed "s/_veloci_/$vel/" < ${TEMPLATE} > ${TARGET}
    for F in gifs_config.ini model.top qm.ndx nve.mdp
    do
        ln -s ${base}/$F .
    done
}

function recordit(){
    k=$1
    surf=$(tail -n1 ${states} | cut -d' ' -f1)
    #dir=$(awk '/1QM/ {if ($7 > 0){ print "T" }else{ print "R" } }' model-test.gro)
    dir=$(if (( $(echo $(grep '1QM' < model-test.gro | cut -c45-53) ">0" | bc -l ) )) ; then echo T ; else echo R; fi)
    echo $k $surf $dir
}

function sample(){
    k=$1
    d=$(mktemp -d)
    (
        cd $d
        buildit $k
        for i in $(seq $N)
        do
            echo "$k: $i/$N"
            runit &> /dev/null
            recordit $k >> ${output}    
        done
    )
    rm -r ${d}
}

#rm ${output} || true

for k in $(seq -w 10.0 ${int} 19.0)
do
   sample $k &
done

wait
