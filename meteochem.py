#!/bin/env python3

import random
import time

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


molecules = ('H2', 'He', 'Ne', 'Ar',
             'CH4', 'NH3', 'H2O', 'H2S',
             'N2',
             'CO', 'CO2', 'O2', 'MgO', 'Al2O3', 'SiO2', 'SO2', 'CaO', 'TiO2', 'Cr2O3', 'MnO',
             'Fe3O4', 'FeS', 'Co3O4', 'CoS', 'Ni2O3', 'NiS', 'Cu', 'Cu2O', 'Cu2S', 'ZnO', 'ZnS',
             
             'HF', 'H3PO4', 'H2SO4', 'HCl',
             'NaOH', 'KOH')

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


mollist = sorted(((mol, moldata(mol)) for mol in molecules), key=lambda x: x[1]['a'], reverse=True)

#for mol, data in mollist:
#    print(f'{mol} => {data!r}')


needmols = 1000

countpermol = {}
totalmols = 0

def nebula_molecules(elemsperdraw, mintw=0):
    elems = dict((n, 0) for n in enames)
    
    while True:
        for e in random.choices(enames, eweights, k=elemsperdraw):
            elems[e] += 1
                
        for mol, data in mollist:
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

                if data['a'] >= mintw:
                    yield (mol, howmany)


cutoff = random.uniform(5,35)
draw = int(random.uniform(1, 10000))
for mol, howmany in nebula_molecules(draw, cutoff):
    if mol in countpermol:
        countpermol[mol] += howmany
    else:
        countpermol[mol] = howmany

    totalmols += howmany

    if totalmols >= needmols:
        break


weightpermol = dict((k, moldata(k)['a'] * v) for k,v in countpermol.items())
totalweight = sum(weightpermol.values())

for mol, molmany in sorted(countpermol.items(), key=lambda x:weightpermol[x[0]], reverse=True):
    print(f"{mol:>7s}: {weightpermol[mol]*100/totalweight:>7.2f}%  ({molmany})")
    
            
            
