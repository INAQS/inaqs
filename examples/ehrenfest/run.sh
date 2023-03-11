#!/bin/bash
set -euo pipefail

readonly NAME=so2_gromacs

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

#../glycine.gro
grompp -f nve.mdp -c so2 -p so2.top -n qm.ndx -o ${NAME}.tpr -maxwarn 1
mdrun -v -deffnm ${NAME} -inaqs inaqs_config.dat
