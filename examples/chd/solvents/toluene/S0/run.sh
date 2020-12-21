#!/bin/bash
source $QCRT/vars.qchem qmmm_external
source /data/home/vale/projects/GQSH/vars.gqsh

set -euo pipefail

readonly VC=/data/home/vale/projects/GQSH/virtualchemistry.org
readonly GIFS=/data/home/vale/projects/GQSH/gifs

# FIXME: better to determine all of these from input or programatically
## Solute parameters
readonly system="1,3-cyclohexadiene in Toluene at 298.15K"
readonly qc_input="../../SF_13_cyclohexadiene_ES_BHHLYP.inp"
readonly RQM="CHD"
readonly qmmm_active="cyclohexadiene"

## Solvent Prameters
readonly solvent_name="toluene_T298.15"
readonly solvent_itp="$VC/GAFF/ITP/toluene.itp"
readonly solvent_pdb="$VC/GAFF/PDB/LIQ/${solvent_name}.pdb.gz"
readonly RSOL="TOL"
readonly solvent_atoms=$(grep -c MOL $solvent_itp)

## GMX name
readonly NAME=13CHD_toluene

# don't litter back-ups all over the place
export GMX_MAXBACKUP=-1

function prepare(){
    # generate .gro and .top files for the QM region
    $GIFS/util/qc2gifs.py -f $qc_input -o ${NAME} -r $RQM -n $qmmm_active -s "$system" --noincludes

    # convert solvent box pdb to gro with residue types indicated and
    # update topology (using __NSOL as a key for solvent molecules)
    # may clobber [defaults] section of passed top
    $GIFS/util/vc2gifs.py -f <(zcat $solvent_pdb) -o $solvent_name -r $RSOL --ff $solvent_itp -p ${NAME}.top
}

function solvate(){
    # make a box
    editconf -f ${NAME}.gro -o ${NAME}.gro -c -d 2 -bt cubic
    
    # pack with Ethanol in a 2 nm shell around the molecule
    #genbox -cp ${NAME}.gro -cs ethanol_T298.15.gro -o ${NAME}.gro -shell 2.0
    genbox -cp ${NAME}.gro -cs ${solvent_name}.gro -o ${NAME}.gro

    # generate a standard index file (${NAME}.ndx) with names like protein, water, system...
    NSolAtoms=$(make_ndx -f ${NAME}.gro -o ${NAME} <<<q | awk "/$RSOL/ {print \$4}")
    NSOL=$(echo $NSolAtoms/$solvent_atoms | bc)
    
    sed -i "s/__NSOL/$NSOL/" ${NAME}.top
}

function run(){
    ulimit -c unlimited
    date
    
    # process all the files into something we can run with
    grompp -f qmmm.mdp -c ${NAME}.gro -p ${NAME}.top -n ${NAME}.ndx -o ${NAME}.tpr

    # do the run; must specify -nt 1
    time mdrun -nt 1 -v -deffnm ${NAME}
    date
}

action=${1:-run}

case "$action" in
    "prepare" | "solvate" | "run")
        $action
        ;;
    "all")
        prepare
        solvate
        run
        ;;
    *)
        echo "valid commands are {prepare,solvate,run,all}; defaults to run for PBS usage"
esac
    
