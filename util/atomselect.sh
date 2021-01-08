#!/bin/bash
set -euo pipefail

readonly exe=$(basename "$0")

function usage(){
    echo "$exe: Use VMD to select atoms from a structure"
    echo "$exe <STRUCTURE> \"SELECTION\""
    echo "  STRUCTURE: some structure legible to VMD"
    echo "  SELECTION: an atomselect string, e.g.: \"not solvent\""
}

## atomselect structure sleection
function atomselect(){
    key="__atomselecting_key"
    vmd -dispdev none "$1" 2>/dev/null \
        <<< "set mysel [ atomselect top \"$2\" ]
             puts \" $key [ \$mysel list ] \"
             quit" |\
                 awk "/$key/ { \$1=\"\"; print \$0 }"
}

#trap usage EXIT

structure=${1:-}

if [ ! -f "$structure" ]
then
    echo "$structure doesn't seem to be a file!"
    usage
    exit
fi

shift
selection="$*"

atomselect "$structure" "$selection" || usage
