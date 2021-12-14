#!/usr/bin/env python3

import fileinput
import numpy as np
import re
import os
from collections import defaultdict
import argparse
import sys
import qc2gifs as qc2g
from pathlib import Path

# convert solvent box pdb to gro with residue types indicated and
# update topology (using __NSOL as a key for solvent molecules)
# may clobber [defaults] section of passed top
# $GIFS/utils/vc2gifs.py -f <(zcat $solvent_pdb) -o $solvent_name -r $RSOL --ff $solvent_itp -p ${NAME}.top


def getArgs():
    parser = argparse.ArgumentParser(description='Convert a virtualchemistry.org solvent box for use with GIFS')
    parser.add_argument("-f", metavar="PDB", required=True,
                        type=Path, help="input solvent box")
    parser.add_argument("--ff", metavar="ITP", required=True,
                        type=Path, help="include forcefield")
    parser.add_argument("-o", metavar="output", required=True,
                        help="name of output .gro/.top files")
    parser.add_argument("-p", metavar="TOP", required=False,
                        type=Path, help="optionally pass a topology file to update")
    parser.add_argument("-r", metavar="residue", default="SOL",
                        help="name for residue; defaults to SOL")
    parser.add_argument("-n", metavar="name", default="VC-SOL",
                        help="system name; defaults to VC-SOL")

    return parser.parse_args()


## PDB format documentation
# http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
# 13 - 16        Atom          name         Atom name.
# 17             Character     altLoc       Alternate location indicator.
# 18 - 20        Residue name  resName      Residue name.
# 22             Character     chainID      Chain identifier.
# 23 - 26        Integer       resSeq       Residue sequence number.
# 27             AChar         iCode        Code for insertion of residues.
# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.


def main():
    args = getArgs()

    itp = readlines(args.ff)
    atypes = getAtoms(itp)

    if args.p is None:
        top = itp
        tfile = args.o + ".out"
    else:
        top = mergeTop(readlines(args.p), itp, args.n, args.r)
        tfile = args.p

    if args.f.suffix != '.gro':
        pdb = readlines(args.f)
        gro = pdb2gro(pdb, atypes, getBox(pdb), args.n, args.r)
    else:
        gro = readlines(args.f)

    gfile = args.o + ".gro"
    with open(gfile, 'w') as out:
        print("\n".join(gro), file=out)
    print(f"wrote {gfile}")

    with open(tfile, "w") as out:
        print("\n".join(top), file=out)
    print(f"wrote {tfile}")

# update topology (using __NSOL as a key for solvent molecules)
# may clobber [defaults] section of passed top


def mergeTop(old_top, itp, solname, residue):
    # take defaults from itp
    # merge atomtypes from itp
    # add moleculetype atoms, bonds, pairs, angles, dihedrals (moleculetype : end)

    top = []
    top += getBlock('defaults', itp)

    top += mergeATypes(itp, old_top)

    # what about the [ nonbond_params ] block?

    top += getBlock('moleculetype', old_top)
    top += getBlock('atoms', old_top)

    top += rename(getBlock('moleculetype', itp), solname)

    # FIXME: need to rename residue
    iatoms = getBlock('atoms', itp)
    for i, line in enumerate(iatoms):
        fields = line.split()
        if len(fields) > 4 and fields[0][0] != ';':
            # insert residue
            fields[3] = residue

            # strip trailing digits from symbol
            # sym = ""
            # for c in fields[4]:
            #     if c.isdigit():
            #         break
            #     else:
            #         sym += c
            # fields[4] = sym
        iatoms[i] = "   ".join(fields)

    top += iatoms
    top += getBlock('bonds', itp)
    top += getBlock('pairs', itp)
    top += getBlock('angles', itp)
    top += getBlock('dihedrals', itp)

    top += getBlock('system', old_top)
    top += getBlock('molecules', old_top)
    top.append(f"{solname} __NSOL")

    return top


def rename(block, name):
    new_block = []
    namere = re.compile(r"[\w-]+\s+\d")
    for line in block:
        if namere.match(line):
            mut = line.split()
            mut[0] = name
            new_block.append(' '.join(mut))
        else:
            new_block.append(line)
    return new_block



def mergeATypes(itp, top):
    atypesblock = ['[ atomtypes ]']

    itpa = getBlock('atomtypes', itp)
    topa = getBlock('atomtypes', top)
    ifloat = r"[e\d\-\.\+]+"
    key = re.compile(f"\w+\s+\w+\s+(\d+\s+)?{ifloat}\s+{ifloat}\s+[ASVD]\s+{ifloat}\s+{ifloat}")

    # FIXME: assumes no different isotopes!
    masses = {}
    for line in getBlock('atoms', itp):
        fields = line.split()
        if len(fields) > 1 and fields[0].isdecimal():
            masses[fields[1]] = float(fields[7])

    for l in itpa + topa:
        if not key.match(l):
            continue
        fields = l.split()
        if not fields[2].isdecimal():
            fields.insert(2, str(guessz(masses[fields[0]])))
        atypesblock.append("    ".join(fields))

    atypesblock.append("\n")

    return atypesblock


# guess the atomic number from the mass
def guessz(m):
    k, diff = None, 1e5
    for mz, z in m2z.items():
        diff_new = abs(m-mz)
        if diff_new < diff:
            diff = diff_new
            k = z
    return k


def pdb2gro(pdb, atypes, box, system, residue):
    gro = []
    atom = re.compile('^ATOM  ')

    N = len(atypes)
    i = 0
    gro.append(system)
    # need to insert total number later

    # freq = defaultdict(lambda : 0)
    for line in pdb:
        if atom.match(line):
            sym = line[13:17].strip()
            # freq[sym] += 1
            # sym += str(freq[sym])
            x = np.array(list(map(float, line[31:55].split())))
            x /= 10  # \AA -> nm
            gro.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" %
                       (int(i/N), residue, sym, i, x[0], x[1], x[2]))
            i += 1
        else:
            continue

    gro.insert(1, str(i))
    gro.append(f"{box[0]} {box[1]} {box[2]}")
    return gro


def getBox(pdb):
    box = re.compile('^CRYST1')
    for l in pdb:
        if box.match(l):
            return np.array(list(map(float, l[7:34].split())))/10


def getAtoms(itp):
    raw = getBlock('atoms', itp)
    key = re.compile(r'\s*\d+.*')

    atypes = {}
    for r in raw:
        if key.match(r):
            n, s = r.split()[0:2]
            atypes[int(n)] = s
    return atypes


def getBlock(block_name, f):
    block = []
    start = re.compile('^\[ *' + block_name + ' *\].*', re.IGNORECASE|re.MULTILINE)
    stop = re.compile('^[\[#].*', re.MULTILINE)
    copy = False
    for line in f:
        if copy and stop.search(line) is not None:
            copy = False
            break
        if start.search(line) is not None:
            copy = True

        if copy:
            block.append(line.strip())

    return block


def readlines(f):
    stripped = []
    with open(f) as lines:
        for l in lines:
            stripped.append(l.strip())
    return stripped


m2z = { 1.00797: 1, 4.00260: 2, 6.941: 3, 9.01218: 4, 10.81: 5,
        12.011: 6, 14.0067: 7, 15.9994: 8, 18.998403: 9, 20.179: 10,
        22.98977: 11, 24.305: 12, 26.98154: 13, 28.0855: 14, 30.97376: 15,
        32.06: 16, 35.453: 17, 39.0983: 19, 39.948: 18, 40.08: 20,
        44.9559: 21, 47.90: 22, 50.9415: 23, 51.996: 24, 54.9380: 25,
        55.847: 26, 58.70: 28, 58.9332: 27, 63.546: 29, 65.38: 30, 69.72:
        31, 72.59: 32, 74.9216: 33, 78.96: 34, 79.904: 35, 83.80: 36,
        85.4678: 37, 87.62: 38, 88.9059: 39, 91.22: 40, 92.9064: 41,
        95.94: 42, 98: 43, 101.07: 44, 102.9055: 45, 106.4: 46, 107.868:
        47, 112.41: 48, 114.82: 49, 118.69: 50, 121.75: 51, 126.9045: 53,
        127.60: 52, 131.30: 54, 132.9054: 55, 137.33: 56, 138.9055: 57,
        140.12: 58, 140.9077: 59, 144.24: 60, 145: 61, 150.4: 62, 151.96:
        63, 157.25: 64, 158.9254: 65, 162.50: 66, 164.9304: 67, 167.26:
        68, 168.9342: 69, 173.04: 70, 174.967: 71, 178.49: 72, 180.9479:
        73, 183.85: 74, 186.207: 75, 190.2: 76, 192.22: 77, 195.09: 78,
        196.9665: 79, 200.59: 80, 204.37: 81, 207.2: 82, 208.9804: 83,
        209: 84, 210: 85, 222: 86, 223: 87, 226.0254: 88, 227.0278: 89,
        231.0359: 91, 232.0381: 90, 237.0482: 93, 238.029: 92, 242: 94,
        243: 95, 247: 97, 247: 96, 250: 102, 251: 98, 252: 99, 255: 108,
        256: 109, 257: 100, 258: 101, 260: 103, 261: 104, 262: 107, 262:
        105, 263: 106, 269: 110, 272: 111, 277: 112}


if __name__ == "__main__":
    main()
