#!/bin/bash
set -euo pipefail

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

readonly NAME=gly_nvt

cp ../gifs_config.ini .

grompp -f nvt.mdp -c nvt.gro -p ../glycine.top -n ../gly.ndx -o ${NAME}.tpr
mdrun -v -deffnm ${NAME}

