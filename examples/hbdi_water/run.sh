#!/bin/bash
set -euo pipefail

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

readonly NAME=HBDI_water

readonly structure=conf.gro
readonly topology=qmmm.top
readonly index=index.ndx
readonly grompp=grompp-slim.mdp
# The index file defines the following groups:
## System, �7, HBD, XXX, NA, Water, SOL, non-Water, Ion, HBD, XXX, NA, SOL_ions, HBD_XXX, System_&_!SOL_ions, a_1-8, QMatoms, System_&_!QMatoms, MMatoms

# process all the files into something we can run with
grompp -f ${grompp} -c ${structure} -p ${topology} -n ${index} -o ${NAME}.tpr

# do the run; must specify -nt 1
mdrun -nt 1 -v -deffnm ${NAME}

