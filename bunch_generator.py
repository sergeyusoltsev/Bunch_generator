#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

slaboffsetz = 0 #slab offset in Angstrom
ss = 10 # slab size x, y, number of point-spheres
ssz = 10 # slab depth, number of point-spheres
pss = 30 # point size (sphere), Angstrom
maxrad = 100
oxarr = []
harr = []

contents = ''
linesflat = ''

with open('NZ64.xyz', 'r') as fh:
    lines = fh.readlines()
    for l in lines[2:]:
        linesflat += l
    # contents += linesflat

mol = np.array(np.genfromtxt('NZ64.xyz', skip_header=2, unpack=True, 
                             dtype=[('name','U5'),('x','<f8'),('y','<f8'),('z','<f8')]))
atoms = []
atnames = []
for i, j in enumerate(mol.T):
    atoms.append([float(j[1]), float(j[2]), float(j[3]), 1])
    atnames.append(str(j[0]))

atoms = np.array(atoms)
atnames = np.array(atnames)

def waterinstance(at = 'He', atom = np.array([[0, 0, 0, 1]]), tr=np.identity(4), raz=0, ray=0, rax=0):
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
    atom = np.dot(atom, Tr)
    return ('%s %.6f %.6f %.6f\n' % 
            (at, atom[0], atom[1], atom[2]))

# ccO = np.array([[0., 0., 0., 1]]) # origin
# ccH1 = np.array([[0.758602, 0., 0.504284, 1]])
# ccH2 = np.array([[0.758602, 0., -0.504284, 1]])
# Translation to starting position in cartesian
tr = np.array([[1,               0,              0,              0],
               [0,               1,              0,              0],
               [0,               0,              1,              0],
               [- ss/2 * pss,    -ss/2 * pss,    -ssz/2 * pss,    1]])

totalatm = 0
fragments = 0

for i in range(ssz):
    for j in range(ss): #y
        for k in range(ss): 
            rax = np.random.rand()*np.pi
            ray = np.random.rand()*np.pi
            raz = np.random.rand()*np.pi#x
            contentsbuff = ''
            totalatmbuff = 0
            for n, atom in enumerate(atoms):
                if np.linalg.norm(np.dot(atom, tr)) < maxrad:
                    contentsbuff += waterinstance(atnames[n][0], atom, tr, raz, ray, rax)
                    totalatmbuff += 1
                else:
                    contentsbuff = r''
                    totalatmbuff = 0
                    break
            contents += contentsbuff
            totalatm += totalatmbuff
            tr[3, 0] += pss
        tr[3, 0] = - ss/2 * pss
        tr[3, 1] += pss
    tr[3, 1] = - ss/2 * pss
    tr[3, 2] += pss
fragments = int(totalatm/len(atnames))

for i in range(1, fragments):
    print('fragment: %i, %i - %i' % (i, totalatm/fragments*(i-1)+1, totalatm/fragments*i))

contents = str(totalatm) + '\n\n' + contents
with open('slab.xyz', 'w') as fh:
    fh.write(contents)
    
print('water slab generated successfully')
