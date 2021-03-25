#!/bin/env python3

import math
import random

## Abundance is relative to silicon=1000

elements = {
    'H':  {'abundance': 27900000,
           'z': 1,
           'a': 1.000},
    'He': {'abundance': 2720000,
           'z': 2,
           'a': 4.003},
    'C':  {'abundance': 10100,
           'z': 6,
           'a': 12.01},
    'N':  {'abundance': 3130,
           'z': 7,
           'a': 14.01},
    'O':  {'abundance': 23800,
           'z': 8,
           'a': 16.00},
    'F':  {'abundance': 1,
           'z': 9,
           'a': 19.00},
    'Ne': {'abundance': 3440,
           'z': 10,
           'a': 20.18},
    'Na': {'abundance': 57,
           'z': 11,
           'a': 2.99},
    'Mg': {'abundance': 1074,
           'z': 12,
           'a': 24.31},
    'Al': {'abundance': 85,
           'z': 13,
           'a': 26.98},
    'Si': {'abundance': 1000,
           'z': 14,
           'a': 28.09},
    'P':  {'abundance': 10,
           'z': 15,
           'a': 30.97},
    'S':  {'abundance': 515,
           'z': 16,
           'a': 32.06},
    'Cl': {'abundance': 5,
           'z': 17,
           'a': 35.45},
    'Ar': {'abundance': 101,
           'z': 18,
           'a': 39.95},
    'K':  {'abundance': 4,
           'z': 19,
           'a': 39.10},
    'Ca': {'abundance': 61,
           'z': 20,
           'a': 40.08},
    'Ti': {'abundance': 2,
           'z': 22,
           'a': 47.87},
    'Cr': {'abundance': 14,
           'z': 24,
           'a': 52.00},
    'Mn': {'abundance': 10,
           'z': 25,
           'a': 54.94},
    'Fe': {'abundance': 900,
           'z': 26,
           'a': 55.85},
    'Co': {'abundance': 2,
           'z': 27,
           'a': 58.93},
    'Ni': {'abundance': 49,
           'z': 28,
           'a': 58.69},
    'Cu': {'abundance': 1,
           'z': 29,
           'a': 63.55},
    'Zn': {'abundance': 1,
           'z': 30,
           'a': 65.38}
    }

molecules = {'H2': {'T3':  13.80},
             'N2': {'T3':  63.15},
             'O2': {'T3':  54.36},
             'F2': {'T3':  53.48},

             'He': {'T3':   2.18},
             'Ne': {'T3':  24.56},
             'Ar': {'T3':  83.81},
             
             'C':  {'T3':  3915},
             'Al': {'T3':  933.47},
             'Mg': {'T3':  923},
             'Si': {'T3':  1687},
             'Ti': {'T3':  1941},
             'Cr': {'T3':  2180},
             'Mn': {'T3':  1519},
             'Fe': {'T3':  1811},
             'Co': {'T3':  1768},
             'Ni': {'T3':  1728},
             'Cu': {'T3':  1358},
             'Zn': {'T3':   693},

             'CH4': {'T3':  90.67},
             'NH3': {'T3': 195.4},
             'H2O': {'T3': 273.16},
             'HF':  {'T3': 190},
             'H2S': {'T3': 187.66},
             'HCl': {'T3': 161.15},

             'CO':    {'T3':  67.9},
             'CO2':   {'T3': 216.58},
             'P2O5':  {'T3': 613},
             'SO2':   {'T3': 197.64},

             'H3PO4': {'T3': 314},

             'MgO':   {'T3': 3125},
             'Al2O3': {'T3': 2345},
             'SiO2':  {'T3': 3220},
             'CaO':   {'T3': 2886},
             'TiO2':  {'T3': 2116},
             'Cr2O3': {'T3': 2708},
             'MnO':   {'T3': 2218},
             'Fe3O4': {'T3': 1870},
             'Co3O4': {'T3': 1168},
             'NiO':   {'T3': 2228},
             'Cu2O':  {'T3': 1505},
             'ZnO':   {'T3': 2247},
             
             'MnS':   {'T3': 1983},
             'FeS':   {'T3': 1467},
             'CoS':   {'T3': 1468},
             'NiS':   {'T3': 1070},
             'Cu2S':  {'T3': 1400},
             'ZnS':   {'T3': 2120},

             'NaF':   {'T3': 1266},
             'NaCl':  {'T3': 1073},
             'Mg2F':  {'T3': 1536},
             'Mg2Cl': {'T3':  987},
             'CaF2':  {'T3': 1691},
             'CaCl2': {'T3': 1046},
             'KF':    {'T3': 1131},
             'KCl':   {'T3': 1040}
             }



def moldata(mol):
    components = {}
    a = 0
    z = 0
    
    s = mol
    while s:
        atom = None
        if s[:2] in elements:
            atom = s[:2]
            s = s[2:]
        elif s[:1] in elements:
            atom = s[:1]
            s = s[1:]
        else:
            raise ValueError(f"Cannot parse molecule {mol}")

        prep = s.lstrip("0123456789");
        howmany = s[:len(s)-len(prep)]
        s = prep
        if howmany:
            howmany = int(howmany)
        else:
            howmany = 1

        if atom in components:
            components[atom] = components[atom] + howmany
        else:
            components[atom] = howmany

        a += elements[atom]['a'] * howmany
        z += elements[atom]['z'] * howmany

    return {'components': components,
            'a': a,
            'z': z}

enames = ()
eweights = ()
for name, data in elements.items():
    print(f'Validating {name} ...')
    enames += (name,)
    eweights += (data['abundance'],)

for name, data in molecules.items():
    print(f'Constructing data sheet for {name} ...')
    data.update(moldata(name))

mollist = sorted(molecules.keys(), key=lambda k: molecules[k]['a'], reverse=True)

for mol in mollist:
    print(f'{mol} => {molecules[mol]!r}')


needmols = random.uniform(1000,20000)

countpermol = {}
totalmols = 0

def nebula_molecules(elemsperdraw, mintw=0, mint3=0):
    elems = dict((n, 0) for n in enames)
    
    while True:
        for e in random.choices(enames, eweights, k=elemsperdraw):
            elems[e] += 1

        #print(repr(elems))
        for mol in mollist:
            data = molecules[mol]
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
                    
cutoff = random.uniform(7,35)
au = random.uniform(2.06,3.27)
#cutoff = random.uniform(0, 35)
#cutoff = 0
#au = random.uniform(0.001,20)
tempatau = 250
temp = math.sqrt(1 / au * tempatau**2)

draw = int(random.uniform(1000, 10000))
for mol, howmany in nebula_molecules(draw, cutoff, temp):
    if mol in countpermol:
        countpermol[mol] += howmany
    else:
        countpermol[mol] = howmany

    totalmols += howmany

    if totalmols >= needmols:
        break


weightpermol = dict((k, moldata(k)['a'] * v) for k,v in countpermol.items())
totalweight = sum(weightpermol.values())

print(f"==== Dist ({au:.2f} AU) ==== Temp ({temp:.1f}K) ==== Flux ({cutoff:.1f} Da) ====")
for mol, molmany in sorted(countpermol.items(), key=lambda x:weightpermol[x[0]], reverse=True):
    print(f"  {mol:>7s}: {weightpermol[mol]*100/totalweight:>7.2f}%  ({molmany})")
    
            
            
