#!/bin/bash
set -euo pipefail

readonly NAME=gly_nve

grompp -f nve.mdp -c ../glycine.gro -p ../glycine.top -n ../gly.ndx -o ${NAME}.tpr
mdrun -v -deffnm ${NAME}
