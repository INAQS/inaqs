import fileinput
import numpy as np
import re
import os
from collections import defaultdict
import argparse
import sys



class Atom:
    S = ''
    x = np.array([])

    def __init__(self, string):
        inp = string.split()
        self.S = inp[0]
        self.x = np.array(list(map(float, inp[1:])))

    def __repr__(self):
        return self.S + " " + " ".join(map(str, self.x))

    def scale(self, a):
        self.x *= a


def getMol(qcfile):
    m = {'molecule': re.compile('^[$]molecule', re.IGNORECASE),
         'end':      re.compile('^[$]end', re.IGNORECASE),
         'bohr':     re.compile('input_bohr.*true', re.IGNORECASE),
         }

    bohr = False
    molecule = []

    # find start of the molecule section
    for line in qcfile:
        if m['bohr'].match(line):
            bohr = True
        if m['molecule'].match(line):
            break

    next(qcfile)  # skip over comment line

    # parse in molecules
    for line in qcfile:
        if m['end'].match(line):
            break

        molecule.append(Atom(line))

    # keep looking for a unit conversion factor
    for line in qcfile:
        if m['bohr'].match(line):
            bohr = True

    scale = 0.10  # \AA -> nm
    if bohr:
        scale *= 0.529177210903  # a_0 -> \AA

    for a in molecule:
        a.scale(scale)

    return molecule


def writeGRO(atoms, name="GIFS", residue="QM", output=sys.stdout):
    """gro file format:

    title string (free format string, optional time in ps after ‘t=’)
    number of atoms (free format integer)
    one line for each atom (fixed format, see below)
    box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z)

    [atoms]
    residue number (5 positions, integer)
    residue name (5 characters)
    atom name (5 characters)
    atom number (5 positions, integer)
    position (in nm, x y z in 3 columns: 8 positions with 3 decimal places)
    velocity (in nm/ps, x y z in 3 columns: 8 positions with 4 decimal places)
    fmt="%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"

    N.B.: We do not inclued velocities
    """
    
    print(name, file=output)
    print(len(atoms), file=output)

    for (i, a) in enumerate(atoms):
        print("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" %
              (1, residue, a.S, i, a.x[0], a.x[1], a.x[2]), file=output)

    print("1 1 1", file=output)


def writeTOP(atoms, name="GIFS", residue="QM", output=sys.stdout):
    # FIXME: want to be able to choose the water topology
    print(f"""
; Include forcefield parameters
#include "oplsaa.ff/forcefield.itp"

; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include water topology; this is for tip3p, but one could use whatever
#include "oplsaa.ff/tip3p.itp"

; Include topology for ions
#include "oplsaa.ff/ions.itp"

[ moleculetype ]
; Name            nrexcl
{name}_QM        3

[ atoms ]
;  nr type         resnr residue atom      cgnr   charge       mass""", file=output)
    for (i, a) in enumerate(atoms, start=1):
        print(fmtOPLSTop(i, a, residue), file=output)

    print(f"""
[ system ]
; Name
{residue}

[ molecules ]
; Compound        #mols
{name}_QM        1
""", file=output)


def fmtOPLSTop(nr, atom, residue, resnr=1):
    #  nr  type       resnr residue  atom    cgnr     charge   mass
    #  1   opls_287      1    GLY      N      1
    #  2   opls_290      1    GLY     H1      1
    # ...
    # or:
    #  1   opls_287      1    GLY      N      1       -0.3    14.0067
    #  2   opls_290      1    GLY     H1      1       0.33      1.008
    # ...

    # Notes: Charge and mass are optional. We do not perturb mass, but
    # set the charge to 0. Charge groups, cgnr, are constructed such
    # that there are no more than 6 atoms in a single charge group (Gromacs' max is 32).
    # The name in the atom field must match that defined in type
    sym = oplsDict[atom.S]
    cgnr = int(nr/(6+1))
    chg = 0.0
    return f"{nr:5} {sym:12} {resnr:5} {residue:7} {atom.S:5} {cgnr:5} {chg:8.3}"


oplsDict = defaultdict(lambda: "opls_???")


def buildOPLSDict():
    gmxdata = os.environ['GMXDATA']

    opls = re.compile("opls_")

    with open(gmxdata + "/top/oplsaa.ff/ffnonbonded.itp") as f:
        for line in f:
            field = line.split()
            if opls.match(field[0]):
                oplsDict[field[1]] = field[0]


def main():
    args = getArgs()
    buildOPLSDict()

    print("reading", args.f)
    with open(args.f) as f:
        molecule = getMol(f)

    with open(args.o + ".gro", 'w') as f:
        writeGRO(molecule, name=args.n, residue=args.r, output=f)
    print("wrote", args.o + ".gro")

    with open(args.o + ".top", 'w') as f:
        writeTOP(molecule, name=args.n, residue=args.r, output=f)
    print("wrote", args.o + ".top")

    print("Be sure to check your topology file for missing atoms")
    
        

def getArgs():
    parser = argparse.ArgumentParser(description='Convert a Q-Chem input file for use with Gromacs under GIFS')
    parser.add_argument("-f", metavar="input", required=True,
                        help="name of Q-Chem input file")
    parser.add_argument("-o", metavar="output", required=True,
                        help="name of output .gro/.top files")
    parser.add_argument("-r", metavar="residue", default="QM",
                        help="name for residue; defaults to QM")
    parser.add_argument("-n", metavar="name", default="GIFS",
                        help="system name; defaults to GIFS")
    return parser.parse_args()


if __name__ == "__main__":
    main()
    # buildOPLSDict()
    # print(oplsDict)
