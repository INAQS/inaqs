#!/bin/bash
set -euo pipefail

readonly NAME=gly_nvt

grompp -f nvt.mdp -c nvt.gro -p ../glycine.top -n ../gly.ndx -o ${NAME}.tpr
mdrun -v -deffnm ${NAME}

