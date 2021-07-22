
#!/usr/bin/env python3

import fileinput
import numpy as np
import re
import os
from collections import defaultdict
import argparse
import sys
from pathlib import Path


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

    print("0 0 0", file=output)


def get_block(block_name, f):
    block = []
    start = re.compile('^\[ *' + block_name + ' *\].*', re.IGNORECASE|re.MULTILINE)
    stop = re.compile('^[\[#].*', re.MULTILINE)
    copy = False
    with open(f) as lines:
        for line in lines:
            line = line.strip()
            if copy and stop.search(line) is not None:
                copy = False
                break
            if start.search(line) is not None:
                copy = True

            if copy:
                block.append(line.strip())
    
    return block


def get_atomtypes(atoms, f):
    ats = set([ ffDict[s2z[atom.S]] for atom in atoms])
    output = []
    with open(f) as lines:
        for line in lines:
            if ats & set(line.split()):
                output.append(line.strip())
    output.append('\n')
    return output


def writeTOP(atoms, name="GIFS", system="GIFS", residue="QM", forcefield=None, noincludes=False, output=sys.stdout):
    # FIXME: want to be able to choose the water topology

    top = []
    
    if noincludes:
        top += get_block('defaults', ffDict['path'] / "forcefield.itp")
        
        top.append('\n[ atomtypes ]')
        top.append('; name bond_type z mass charge ptype sigma espilon')
        top += get_atomtypes(atoms, ffDict['path'] / "ffnonbonded.itp")

    else:
        top.append('; Include forcefield parameters')
        top.append(f'#include "{forcefield}/forcefield.itp"')
    
    top.append("\n[ moleculetype ]")
    top.append("; Name         nrexcl")
    top.append(f"{name}        3")

    top.append("\n[ atoms ]")
    top.append(";  nr type         resnr residue atom      cgnr   charge       mass")

    for (i, a) in enumerate(atoms, start=1):
        top.append(fmtTopAtoms(i, a, residue))

    top.append("\n[ system ]")
    top.append("; Name")
    top.append(f"QMMM {system}")
    top.append("\n[ molecules ]")
    top.append("; Compound        #mols")
    top.append(f"{name}        1")

    print("\n".join(top), file=output)



def fmtTopAtoms(nr, atom, residue, resnr=1):
    #  nr  type       resnr residue  atom    cgnr     charge   mass
    #  1   opls_287      1    GLY      N      1
    #  2   opls_290      1    GLY     H1      1
    # ...
    # or:
    #  1   opls_287      1    GLY      N      1       -0.3    14.0067
    #  2   opls_290      1    GLY     H1      1       0.33      1.008
    # ...

    # Notes: Charge and mass are optional. We do not perturb mass, but
    # set the charge to 0. The name in the atom field must match that
    # defined in type. Charge groups, cgnr, *were previously*
    # constructed such that there are no more than 6 atoms in a single
    # charge group (Gromacs' max is 32). However, this produced errors
    # (blowing-up) with some solvents. Following qforce [1], we simply
    # place each atom in its own charge group.

    # [1]: https://github.com/selimsami/qforce/blob/master/forcefield.py#L168

    sym = ffDict[s2z[atom.S]]
    cgnr = nr
    chg = 0.0
    return f"{nr:5} {sym:12} {resnr:5} {residue:7} {atom.S:5} {cgnr:5} {chg:8.3}"


def genffD(forcefield):

    ff = Path(forcefield)
    if not ff.exists():
        ff = Path(os.environ['GMXDATA']) / "top" / (forcefield)

    ffDict['path']=ff

    with open (ff / "ffnonbonded.itp") as f:
        o = 0  # special case offset for OPLS-AA
        if forcefield == "oplsaa.ff":
            o = 1
        for line in f:
            field = line.split()
            if len(field) > 4 + o and field[4+o] == 'A':
                ffDict[int(field[1+o])] = field[0]
                    

def main():
    args = getArgs()
    genffD(args.ff)

    print("reading", args.f)
    with open(args.f) as f:
        molecule = getMol(f)

    with open(args.o + ".gro", 'w') as f:
        writeGRO(molecule, name=args.n, residue=args.r, output=f)
    print("wrote", args.o + ".gro")

    with open(args.o + ".top", 'w') as f:
        writeTOP(molecule, name=args.n, system=args.s, residue=args.r, forcefield=args.ff, noincludes=args.noincludes, output=f)
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
                        help="molecule name; defaults to GIFS")
    parser.add_argument("-s", metavar="name", default="GIFS",
                        help="system name; defaults to GIFS")
    parser.add_argument("--ff", metavar="forcefield", default="oplsaa.ff",
                        help="forcefield to use; defaults to oplsaa.ff")
    parser.add_argument("--noincludes", default=False, action='store_true',
                        help="generate self-contained topology with no includes")
    return parser.parse_args()


ffDict = defaultdict(lambda: print("Warning: unknown atom in forcefield", file=sys.stderr) or "ff_???")

s2z = {
    'H':   1, 'He':  2, 'Li':  3, 'Be':  4, 'B':   5, 'C':   6,
    'N':   7, 'O':   8, 'F':   9, 'Ne': 10, 'Na': 11, 'Mg': 12,
    'Al': 13, 'Si': 14, 'P':  15, 'S':  16, 'Cl': 17, 'Ar': 18,
    'K':  19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V':  23, 'Cr': 24,
    'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y':  39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
    'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
    'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I':  53, 'Xe': 54,
    'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
    'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
    'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72,
    'Ta': 73, 'W':  74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
    'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84,
    'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
    'Pa': 91, 'U':  92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96,
    'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102,
    'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
    'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114,
    'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118,
}

z2s = {v:k for k,v in s2z.items()}

if __name__ == "__main__":
    main()
    # genffD("oplsaa.ff")
    # buildffDict()
    # print(ffDict)
