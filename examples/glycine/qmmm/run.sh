#!/bin/bash
set -euo pipefail

readonly NAME=gly_qmmm

export GMX_MAXBACKUP=-1

cp ../glycine.top .
editconf -f qmmm.gro -o qmmm-solvated.gro -c -d 2 -bt cubic
genbox -cp qmmm-solvated.gro -cs spc216.gro -o qmmm-solvated.gro -shell 1.0 -p glycine.top
make_ndx -f qmmm-solvated.gro -o ${NAME} <<<q

grompp -f qmmm.mdp -c qmmm-solvated.gro -p glycine.top -n ${NAME}.ndx -o ${NAME}.tpr
#mdrun -v -deffnm ${NAME}

