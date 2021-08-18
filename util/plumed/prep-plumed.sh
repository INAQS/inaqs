#!/bin/bash
set -euo pipefail

PK=${PLUMED_KERNEL:-}

if [ -z "$PK" ]
then
    echo "Plumed seems not to be installed; PLUMED_KERNEL is not set"
    exit 1
fi

HERE=$(readlink -f "$(dirname "$0")")

PLUMED_LIB=$(readlink -f "$(dirname "$PLUMED_KERNEL")")

echo "$PLUMED_LIB/patches"
echo "$PLUMED_LIB/plumed/patches"

if [ -d "$PLUMED_LIB/plumed/patches" ]
then
    PLUMED_PATCHES="$PLUMED_LIB/plumed/patches"
elif [ -d "$PLUMED_LIB/../../patches" ]
then
    PLUMED_PATCHES="$PLUMED_LIB/patches"
else
    echo "can't find patches directory for plumed install"
    exit 1
fi
     
GMX=gromacs-4.6.5

cp -at "$PLUMED_PATCHES" "$HERE"/$GMX.{config,diff}

echo "Patches for $GMX added to your plumed install at"
echo "$PLUMED_PATCHES"
echo "You should be able to go to the gromacs install in gifs, execute"
echo "plumed patch --patch --engine $GMX"
echo "and then recompile gromacs."
