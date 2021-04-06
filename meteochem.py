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

# T3: Triple point temperature (K), where available.  Melting point is used when unavailable.
# Hf: Standard enthalpy of formation (kJ/mol)

molecules = {'H2': {'T3':  13.80, 'Hf': 0.0},
             'N2': {'T3':  63.15, 'Hf': 0.0},
             'O2': {'T3':  54.36, 'Hf': 0.0},
             'F2': {'T3':  53.48, 'Hf': 0.0},

             'He': {'T3':   2.18, 'Hf': 0.0},
             'Ne': {'T3':  24.56, 'Hf': 0.0},
             'Ar': {'T3':  83.81, 'Hf': 0.0},
             
             'C':  {'T3':  3915, 'Hf': 0.0},
             'Al': {'T3':  933.47, 'Hf': 0.0},
             'Mg': {'T3':  923, 'Hf': 0.0},
             'Si': {'T3':  1687, 'Hf': 0.0},
             'Ti': {'T3':  1941, 'Hf': 0.0},
             'Cr': {'T3':  2180, 'Hf': 0.0},
             'Mn': {'T3':  1519, 'Hf': 0.0},
             'Fe': {'T3':  1811, 'Hf': 0.0},
             'Co': {'T3':  1768, 'Hf': 0.0},
             'Ni': {'T3':  1728, 'Hf': 0.0},
             'Cu': {'T3':  1358, 'Hf': 0.0},
             'Zn': {'T3':   693, 'Hf': 0.0},

             'CH4': {'T3':  90.67, 'Hf': -74.6},
             'NH3': {'T3': 195.4,  'Hf': -46},
             'H2O': {'T3': 273.16, 'Hf': -291.83},
             'HF':  {'T3': 190,    'Hf': -13.66},
             'H2S': {'T3': 187.66, 'Hf': -21},
             'HCl': {'T3': 161.15, 'Hf': -92.31},

             'CO':    {'T3':  67.9,  'Hf': -110.5},
             'CO2':   {'T3': 216.58, 'Hf': -393.5},
             'P2O5':  {'T3': 613,    'Hf': -1452},
             'SO2':   {'T3': 197.64, 'Hf': -296.81},

             'H3PO4': {'T3': 314, 'Hf': -1271.7},

             'MgO':   {'T3': 3125, 'Hf': -601.6},
             'Al2O3': {'T3': 2345, 'Hf': -1675.7},
             'SiO2':  {'T3': 3220, 'Hf': -911},
             'CaO':   {'T3': 2886, 'Hf': -635},
             'TiO2':  {'T3': 2116, 'Hf': -945},
             'Cr2O3': {'T3': 2708, 'Hf': -1128},
             'MnO':   {'T3': 2218, 'Hf': -385},
             'Fe3O4': {'T3': 1870, 'Hf': -1120.89},
             'CoO':   {'T3': 2206, 'Hf': -237.74},
             'NiO':   {'T3': 2228, 'Hf': -240},
             'Cu2O':  {'T3': 1505, 'Hf': -170},
             'ZnO':   {'T3': 2247, 'Hf': -350.5},
             
             'MnS':   {'T3': 1983, 'Hf': -345.72},
             'FeS':   {'T3': 1467, 'Hf': -101.67},
             'CoS':   {'T3': 1468, 'Hf': -190},  # Estimated?
             'NiS':   {'T3': 1070, 'Hf': -87.86},
             'Cu2S':  {'T3': 1400, 'Hf': -120},  # Estimated, badly?
             'ZnS':   {'T3': 2120, 'Hf': -204.6},

             'NaF':   {'T3': 1266, 'Hf': -573.6},
             'NaCl':  {'T3': 1073, 'Hf': -411.12},
             'Mg2F':  {'T3': 1536, 'Hf': -1124.2},
             'Mg2Cl': {'T3':  987, 'Hf': -641.1},
             'CaF2':  {'T3': 1691, 'Hf': -1225.91},
             'CaCl2': {'T3': 1046, 'Hf': -795.42},
             'KF':    {'T3': 1131, 'Hf': -568.61},
             'KCl':   {'T3': 1040, 'Hf': -436},

             # Orthosilicates - Olivines
             'Mg2SiO4':  {'T3': 2163, 'Hf': -2056},
             'Fe2SiO4':  {'T3': 1473, 'Hf': -1379},
             }



def moldata(mol):
    components = {}
    a = 0
    n = 0
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
        n += howmany
        z += elements[atom]['z'] * howmany

    return {'components': components,
            'a': a,
            'n': n,
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
    data['aepa'] = data['Hf'] / data['n']
    data['aepb'] = 0.0 if data['n'] == 1 else data['Hf'] / (data['n'] - 1.0)

mollist = sorted(molecules.keys(), key=lambda k: (molecules[k]['aepa'], -molecules[k]['a']))

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


choice = random.randint(0,2)
if choice == 0:
    #cutoff = random.uniform(0,35)
    au = random.uniform(2.06,3.27)
elif choice == 1:
    #cutoff = random.uniform(0, 5)
    au = random.uniform(30, 50)
elif choice == 2:
    au = random.uniform(20000,50000)
    
    
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


weightpermol = dict((k, moldata(k)['a'] * v) for k,v in countpermol.items())
totalweight = sum(weightpermol.values())

print(f"==== Dist ({au:.2f} AU) ==== Temp ({temp:.1f}K) ==== Flux ({cutoff:.1f} Da) ====")
for mol, molmany in sorted(countpermol.items(), key=lambda x:weightpermol[x[0]], reverse=True):
    print(f"  {mol:>7s}: {weightpermol[mol]*100/totalweight:>7.2f}%  ({molmany})")
    
            
            
