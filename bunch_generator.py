#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os, argparse

parser = argparse.ArgumentParser(description='Generate molecule with explicit solvation environment')
parser.add_argument('-mol','--mol', metavar='molecule', type=str, default='molecule.xyz',
                    help='name of molecule file')
parser.add_argument('-sol','--sol', metavar='solvent', type=str, default='solvent.xyz',
                    help='name of solvent file')
parser.add_argument('-ss','--ss', metavar='N', type=int, default=35,
                    help='slab size in Angstroem')
parser.add_argument('-pss','--pss', metavar='N', type=float, default=6.5,
                    help='size of single molecule occupied space in Angstroem')
parser.add_argument('-mr','--mr', metavar='N', type=float, default=60,
                    help='cutoff radius in Angstroem')
parser.add_argument('-cr','--cr', metavar='N', type=float, default=3.5,
                    help='collide solvent-molecule radius in Angstroem')
parser.add_argument('-lr','--lr', metavar='N', type=float, default=15,
                    help='layer radius in Angstroem, set huge number to generate sphere of --mr')
parser.add_argument('-rt','--rt', metavar='type', type=str, default='gen',
                    help='Runtype, use gen for generation or trim for trimming with -nat set to number of atoms in the beginning of .xyz as comparator')
parser.add_argument('-nat','--nat', metavar='N', type=int, default=270,
                    help='number of atoms in the beginning of .xyz to be used as comparator in trimming procedure')
parser.add_argument('-nats','--nats', metavar='N', type=int, default=6,
                    help='number of atoms in each solvent molecule')

args = parser.parse_args()

molname = args.mol
solvname = args.sol
ss = args.ss
pss = args.pss
maxrad = args.mr
colliderad = args.cr
layerrad = args.lr
runtype = args.rt
natoms = args.nat
solvnatoms = args.nats

path = os.path.abspath(__file__)[:-len(os.path.basename(__file__))]

# Starting position in cartesian

tr = np.array([[1,               0,              0,              0],
               [0,               1,              0,              0],
               [0,               0,              1,              0],
               [- ss/2 * pss,    -ss/2 * pss,    -ss/2 * pss,    1]])

glob_contents = '' # this is main global we will write to generated xyz file
glob_fragments = []
glob_atom_amnt = 0

path = os.path.abspath(__file__)[:-len(os.path.basename(__file__))]

# imports moleciles in path with sourcename (optionally according to the slice)

def importxyz(sourcename, path, slicebeg=0, sliceend=None):
    mol = np.array(np.genfromtxt(str(path) + sourcename, skip_header=2, unpack=True, 
                                 dtype=[('name','U5'),('x','<f8'),('y','<f8'),('z','<f8')]))
    atoms = []
    atnames = []
    for i, j in enumerate(mol.T):
        atoms.append([float(j[1]), float(j[2]), float(j[3]), 1])
        atnames.append(str(j[0]))
    atoms = np.array(atoms[slicebeg:sliceend])
    atnames = np.array(atnames[slicebeg:sliceend])
    amnt = len(atnames[slicebeg:sliceend])
    return atoms, atnames, amnt

# generates STRINGS from xyz files containing valid xyz form AND fragment data

def generatefromxyz(atoms, atnames, amnt):
    contents = ''
    frag = ''
    for n, atom in enumerate(atoms):
        contents += instantiate(atnames[n], atom)
    frag = len(atoms)
    return contents, frag

# moves single atom around according to translation vector and angles in radians

def move(atom = np.array([[0, 0, 0, 1]]), tr=np.identity(4), raz=0, ray=0, rax=0):
    Rz = np.array([[np.cos(raz),     np.sin(raz),    0,  0],
                   [-np.sin(raz),    np.cos(raz),    0,  0],
                   [0,               0,              1,  0],
                   [0,               0,              0,  1]])
    Ry = np.array([[-np.sin(ray),    np.cos(ray),    0,  0],
                   [0,               0,              1,  0],
                   [np.cos(ray),     np.sin(ray),    0,  0],
                   [0,               0,              0,  1]])
    Rx = np.array([[0,               0,              1,  0],
                   [np.cos(rax),     np.sin(rax),    0,  0],
                   [-np.sin(rax),    np.cos(rax),    0,  0],
                   [0,               0,              0,  1]])
    Tr = np.dot(Rz, np.dot(Ry, np.dot(Rx, tr)))
    return np.dot(atom, Tr)

# returns text form of atom in xyz format

def instantiate(at = 'He', atom = np.array([[0, 0, 0, 1]])):
    return ('%s %.6f %.6f %.6f\n' % 
            (at, atom[0], atom[1], atom[2]))

# checks if one atom is too close to any of the atoms of the collider

def checkcollision(atom, collider, colliderad):
    for i, j in enumerate(collider):
        distance = np.sqrt((j[0] - atom[0])**2 + (j[1] - atom[1])**2 + (j[2] - atom[2])**2)
        if distance < colliderad:
            return False
    return True

# checks if atom is close enough to any of the atoms of the collider

def checktoofarfrom(atom, collider, layerrad):
    for i, j in enumerate(collider):
        distance = np.sqrt((j[0] - atom[0])**2 + (j[1] - atom[1])**2 + (j[2] - atom[2])**2)
        if distance > layerrad:
            closenoughflag = False
        else:
            closenoughflag = True
            break
    return closenoughflag

# checks if new placement is valid and returns molecule data and fragment data if OK

def generateifvalid(atoms, atnames, frag, collider, tr, rax, ray, raz):
    contents = r''
    frag = frag
    for n, atom in enumerate(atoms):
        current = move(atom, tr, rax, ray, raz)
        if np.linalg.norm(current, ord=2) < maxrad and checkcollision(current, collider, colliderad) and checktoofarfrom(current, collider, layerrad):
            contents += instantiate(atnames[n], current)
        else:
            return False
    return contents, frag

# programs starts here

def runtype_gen(molname, solvname, ss, pss, maxrad, colliderad, layerrad,
                tr, glob_contents, glob_fragments, glob_atom_amnt, path):
    base_atoms, base_atnames, base_amnt = importxyz(molname, path)
    contents, frag = generatefromxyz(base_atoms, base_atnames, base_amnt)
    glob_fragments.append(str(frag))
    glob_atom_amnt += frag
    glob_contents += contents
    
    surr_atoms, surr_atnames, surr_amnt = importxyz(solvname, path)
    print(surr_atnames)
    for i in range(ss):
        for j in range(ss): #y
            for k in range(ss): 
    
                rax = np.random.rand()*2*np.pi
                ray = np.random.rand()*2*np.pi
                raz = np.random.rand()*2*np.pi#x
    
                surr_contents, surr_frag = generatefromxyz(surr_atoms, surr_atnames, surr_amnt)
                try:
                    append_contents, append_frag = generateifvalid(surr_atoms, surr_atnames, surr_frag, base_atoms, tr, rax, ray, raz)
                except:
                    continue
                finally:
                    tr[3, 0] += pss
                glob_fragments.append(str(append_frag))
                glob_atom_amnt += append_frag
                glob_contents += append_contents
    
            tr[3, 0] = - ss/2 * pss
            tr[3, 1] += pss
    
        tr[3, 1] = - ss/2 * pss
        tr[3, 2] += pss
    
    print(glob_fragments)
    
    with open(str(path) + '.xcontrol', 'w') as fh:
        counter = 0
        fh.write('$split\n')
        for n, frag in enumerate(glob_fragments):
            fh.write(' fragment%i:%i-%i\n' % (n + 1, counter + 1, counter + int(frag)))
            counter += int(frag)
    
    glob_contents = str(glob_atom_amnt) + '\n\n' + glob_contents
    with open(str(path) + 'slab.xyz', 'w') as fh:
        fh.write(glob_contents)
    
    print('Molecular bunch of %i atoms generated successfully in slab.xyz, .xcontrol file contains all the fragments' % (glob_atom_amnt))

def runtype_trim(molname, solvname, ss, pss, maxrad, colliderad, layerrad, 
                 natoms, path, tr, glob_contents, glob_fragments, glob_atom_amnt, solvnatoms):
    base_atoms, base_atnames, base_amnt = importxyz(molname, path, sliceend = natoms)
    contents, frag = generatefromxyz(base_atoms, base_atnames, base_amnt)
    glob_fragments.append(str(frag))
    glob_atom_amnt += natoms
    glob_contents += contents

    for i in range(0, 999999):
        leaveme = False
        slicebeg=(natoms+i*solvnatoms)
        sliceend=(natoms+(i+1)*solvnatoms)
        frag_atoms, frag_atnames, frag_amnt = importxyz(molname, path, slicebeg, sliceend)
        if len(frag_atoms)==0:
            break
        collider = base_atoms
        for n, atom in enumerate(frag_atoms):
            if checktoofarfrom(atom, collider, layerrad):
                leaveme = True
            else:
                leaveme = False
                break
        if leaveme == True:
            contents, frag = generatefromxyz(frag_atoms, frag_atnames, frag_amnt)
            glob_fragments.append(str(frag))
            glob_atom_amnt += frag
            glob_contents += contents
    
    print(glob_fragments)
    
    with open(str(path) + '.xcontrol', 'w') as fh:
        counter = 0
        fh.write('$split\n')
        for n, frag in enumerate(glob_fragments):
            fh.write(' fragment%i:%i-%i\n' % (n + 1, counter + 1, counter + int(frag)))
            counter += int(frag)
    
    glob_contents = str(glob_atom_amnt) + '\n\n' + glob_contents
    with open(str(path) + 'slabtrimmed.xyz', 'w') as fh:
        fh.write(glob_contents)
    
    print('Molecular bunch trimmed to %i atoms to slabtrimmed.xyz, .xcontrol file contains all the fragments' % (glob_atom_amnt))
    
    
if runtype == 'gen':
    runtype_gen(molname, solvname, ss, pss, maxrad, colliderad, layerrad,
                tr, glob_contents, glob_fragments, glob_atom_amnt, path)
else:
    runtype_trim(molname, solvname, ss, pss, maxrad, colliderad, layerrad,
                 natoms, path, tr, glob_contents, glob_fragments, glob_atom_amnt,
                 solvnatoms)