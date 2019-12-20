#!/bin/bash
set -euo pipefail

readonly NAME=gly_nve

#../glycine.gro
grompp -f nve.mdp -c nve.gro  -p ../glycine.top -n ../gly.ndx -o ${NAME}.tpr
mdrun -v -deffnm ${NAME}
