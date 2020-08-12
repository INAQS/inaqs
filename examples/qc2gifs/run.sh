#!/bin/bash
set -euo pipefail

readonly NAME=closs1_nve

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

python3 ../../util/qc2gifs.py -o ${NAME} -f ./closs1.in
make_ndx -f ${NAME}.gro -o ${NAME} <<< q

grompp -f nve.mdp -c ${NAME}.gro  -p ${NAME}.top -n ${NAME}.ndx -o ${NAME}.tpr
mdrun -v -deffnm ${NAME}
