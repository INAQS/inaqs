#!/bin/bash
set -euo pipefail

readonly NAME=gly_nve

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

#../glycine.gro
grompp -f nve.mdp -c nve.gro  -p ../glycine.top -n ../gly.ndx -o ${NAME}.tpr
mdrun -v -deffnm ${NAME}
