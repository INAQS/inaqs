#!/bin/bash
set -euo pipefail

readonly NAME=gly_qmmm

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

# pull in gifs_config
[ ! -f gifs_config.ini ] && cp ../gifs_config.ini .

# the topology file will be modified so bring in our own
cp ../glycine.top .

# center the molecule in a large box
editconf -f qmmm.gro -o qmmm-solvated.gro -c -d 2 -bt cubic

# pack with waters in a 1 nm shell around the molecule
genbox -cp qmmm-solvated.gro -cs spc216.gro -o qmmm-solvated.gro -shell 1.0 -p glycine.top

# generate a standard index file with names like protein, water, system...
make_ndx -f qmmm-solvated.gro -o ${NAME} <<<q

# process all the files into something we can run with
grompp -f qmmm.mdp -c qmmm-solvated.gro -p glycine.top -n ${NAME}.ndx -o ${NAME}.tpr

# do the run; must specify -nt 1
mdrun -nt 1 -v -deffnm ${NAME}

