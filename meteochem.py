#!/bin/env python3

import math
import moldata
import random




needmols = random.uniform(1000,20000)

countpermol = {}
totalmols = 0

def nebula_molecules(elemsperdraw, mintw=0, mint3=0):
    elems = dict((n, 0) for n in moldata.enames)
    
    while True:
        for e in random.choices(moldata.enames, moldata.eweights, k=elemsperdraw):
            elems[e] += 1

        #print(repr(elems))
        for mol in moldata.mollist:
            data = moldata.molecules[mol]
            howmany = None
            for e, n in data['components'].items():
                hm = elems[e] // n
                if hm:
                    if howmany is None or howmany > hm:
                        howmany = hm
                else:
                    howmany = 0
                    break

            if howmany:
                for e, n in data['components'].items():
                    elems[e] -= howmany * n

                if data['a'] >= mintw and data['T3'] >= mint3:
                    yield (mol, howmany)


# E <-> T**4    # Stefan-Boltzmann Law  (Is proportional to)
# E <-> 1/d**2   or   d**2 <-> 1/E   or   d <-> sqrt(1/E)
# d <-> sqrt(1/T**4)   or  d <-> 1/T**2
# (d1 * T1**2 == d2 * T2**2)   or   T2**2 == d1 / d2 * T1**2


## choice = random.randint(0,3)
## if choice == 0:
##     au = random.uniform(0.1, 2.0)
## if choice == 1:
##     #cutoff = random.uniform(0,35)
##     au = random.uniform(2.06,3.27)
## elif choice == 2:
##     #cutoff = random.uniform(0, 5)
##     au = random.uniform(30, 50)
## elif choice == 3:
##     au = random.uniform(20000,50000)
sigma = (2.06 + 3.27) / 2.0
left = random.gauss(0,sigma)
right = random.gauss(0,sigma)
au = math.sqrt(left * left + right * right)
    
#cutoff = random.uniform(0, 35)
cutoff = 0
#au = random.uniform(0.001,20)
tempatau = 250
temp = max(math.sqrt(1 / au * tempatau**2), 3.2)

draw = int(random.uniform(1000, 10000))
for mol, howmany in nebula_molecules(draw, cutoff, temp):
    if mol in countpermol:
        countpermol[mol] += howmany
    else:
        countpermol[mol] = howmany

    totalmols += howmany

    if totalmols >= needmols:
        break


weightpermol = dict((k, moldata.molecules[k]['a'] * v) for k,v in countpermol.items())
totalweight = sum(weightpermol.values())

print(f"==== Dist ({au:.2f} AU) ==== Temp ({temp:.1f}K) ==== Flux ({cutoff:.1f} Da) ====")
for mol, molmany in sorted(countpermol.items(), key=lambda x:weightpermol[x[0]], reverse=True):
    print(f"  {mol:>7s}: {weightpermol[mol]*100/totalweight:>7.2f}%  ({molmany})")
    
            
            
