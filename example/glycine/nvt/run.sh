#!/bin/bash
set -euo pipefail

readonly NAME=gly_minimize

grompp -f minim.mdp -c glycine.gro -p ../glycine.top -n ../gly.ndx -o ${NAME}.tpr
mdrun -v -deffnm ${NAME}

echo
echo "Your minimized structure should now be in ${NAME}.gro"
