#!/bin/bash
set -euo pipefail

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

readonly NAME=HBDI_water

readonly structurein=conf.gro
readonly structure=conf-centered.gro
readonly topology=qmmm.top
readonly index=index.ndx
# The index file defines the following groups:
## System, �7, HBD, XXX, NA, Water, SOL, non-Water, Ion, HBD, XXX, NA, SOL_ions, HBD_XXX, System_&_!SOL_ions, a_1-8, QMatoms, System_&_!QMatoms, MMatoms

# center the molecule in a large box
editconf -f ${structurein} -o ${structure} -c

# process all the files into something we can run with
#grompp -f qmmm.mdp -c qmmm-solvated.gro -p glycine.top -n ${NAME}.ndx -o ${NAME}.tpr

# do the run; must specify -nt 1
#mdrun -nt 1 -v -deffnm ${NAME}

