#!/bin/bash
set -euo pipefail

readonly NAME=model-test

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

#../glycine.gro
grompp -f nve.mdp -c model -p model.top -n qm.ndx -o ${NAME}.tpr -maxwarn 1
mdrun -v -deffnm ${NAME}
