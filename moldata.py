#!/bin/env python3

import numpy as np;
import random;

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

             'Na2O':  {'T3': 1405, 'Hf': -416},
             'MgO':   {'T3': 3125, 'Hf': -601.6},
             'Al2O3': {'T3': 2345, 'Hf': -1675.7},
             'SiO2':  {'T3': 3220, 'Hf': -911},
             'K2O':   {'T3': 1010, 'Hf': -363.17},
             'CaO':   {'T3': 2886, 'Hf': -635},
             'TiO2':  {'T3': 2116, 'Hf': -945},
             'Cr2O3': {'T3': 2708, 'Hf': -1128},
             'MnO':   {'T3': 2218, 'Hf': -385},
             'FeO':   {'T3': 1650, 'Hf': -272.04},
             'Fe2O3': {'T3': 1812, 'Hf': -824.20},
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

             'H3PO4': {'T3': 314,    'Hf': -1271.7},
             'H2SO4': {'T3': 283.46, 'Hf': -814},

             'NaOH':   {'T3': 596, 'Hf': -425.8},
             'MgO2H2': {'T3': 623, 'Hf': -924.7},
             'KOH':    {'T3': 633, 'Hf': -425.8},
             'CaO2H2': {'T3': 853, 'Hf': -987},

             'NaF':   {'T3': 1266, 'Hf': -573.6},
             'NaCl':  {'T3': 1073, 'Hf': -411.12},
             'MgF2':  {'T3': 1536, 'Hf': -1124.2},
             'MgCl2': {'T3':  987, 'Hf': -641.1},
             'KF':    {'T3': 1131, 'Hf': -568.61},
             'KCl':   {'T3': 1040, 'Hf': -436},
             'CaF2':  {'T3': 1691, 'Hf': -1225.91},
             'CaCl2': {'T3': 1046, 'Hf': -795.42},

             }



def process(mol):
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
    enames += (name,)
    eweights += (data['abundance'],)

for name, data in molecules.items():
    data.update(process(name))
    data['aepa'] = data['Hf'] / data['n']
    data['aepb'] = 0.0 if data['n'] == 1 else data['Hf'] / (data['n'] - 1.0)

mollist = sorted(molecules.keys(), key=lambda k: (molecules[k]['aepa'], -molecules[k]['a']))

def silica_content(data):
    if 'Si' in data['components'] and 'O' in data['components']:
        m = min(data['components']['Si'], data['components']['O'] // 2)
        w = molecules['SiO2']['a'] * m
        return w / data['a']
    else:
        return 0

def silica_classify(frac):
    return("Non-silicous"        if frac <= 0.0 else
           "Ultramafic"          if frac < 0.45 else
           "Mafic"               if frac < 0.52 else
           "Intermediate"        if frac < 0.63 else
           "Intermediate-felsic" if frac < 0.69 else
           "Felsic")
        

def balance(reagents, products):
    atomsA = dict();
    rpsize = len(reagents) + len(products);

    for index, mol in enumerate(reagents):
        for atom, howmany in molecules[mol]['components'].items():
            if atom not in atomsA:
                atomsA[atom] = [0] * rpsize
            atomsA[atom][index] = howmany

    for index, mol in enumerate(products, start=len(reagents)):
        for atom, howmany in molecules[mol]['components'].items():
            if atom not in atomsA:
                atomsA[atom] = [0] * rpsize
            atomsA[atom][index] = -howmany

    
    print(repr(atomsA));      
    
    # R0 * nr0 + R1 * nr1 = P0 * np0 + P1 * np1
    # R0 * mr0 + R1 * mr1 = P0 * mp0 + P1 * mp1


## In the form, a set indicates alternation, and a tuple indicates serialization
_MO = (frozenset(('MgO', 'MnO', 'FeO', 'CoO', 'NiO', 'ZnO')),)
_SIO2 = ('SiO2',)
_XO = (frozenset(('MgO', 'CaO', 'MnO', 'FeO', 'CoO', 'NiO', 'ZnO')),)
_A2O = (frozenset(('Na2O', 'K2O')),)
_X2O3 = (frozenset(('Al2O3', 'Cr2O3', 'Fe2O3')),)
_AOH = (frozenset(('NaOH','KOH')),)
_AAOHOH = (frozenset(((_AOH + _AOH), 'MgO2H2', 'CaO2H2')),)
_AC = (frozenset(('NaOH', 'NaF', 'NaCl', 'KOH', 'KF', 'KCl')),)
_AACC = (frozenset(((_AC + _AC), 'MgO2H2', 'MgF2', 'MgCl2', 'CaO2H2', 'CaF2', 'CaCl2')),)


minerals = (
    {'name': 'Olivine',
     'form': 2*_MO + _SIO2 },
    {'name': 'Titanite',
     'form': ('CaO', 'TiO2') + _SIO2},
    {'name': 'Plagioclase',
     'form': (frozenset((('CaO', 'CaO', 'Al2O3'), ('Na2O',) + 2*_SIO2)),) + ('Al2O3',) + 4*_SIO2 },
    {'name': 'Pyroxene',
     'form': (frozenset((4*_XO, _A2O + _X2O3)),) +  4*_SIO2 },
    {'name': 'Amphibole',
     'form': _AAOHOH + (frozenset((4*_XO, 2*_A2O + _X2O3)),) + 2 * _XO + 8 * _SIO2 },
    {'name': 'Mica',
     'form': _AACC + _AACC + 2*_XO + 2*_X2O3 + 6 * _SIO2 },
    {'name': 'Orthoclase',
     'form': _A2O + ('Al2O3',) + 6 * _SIO2 },
    {'name': 'Quartz',
     'form': _SIO2 },
);


def drawfrom(mineralform, molmany):
    if isinstance(mineralform, str):
        if mineralform in molmany and molmany[mineralform] > 0:
            molmany[mineralform] -= 1
            return (mineralform,)
        else:
            return None
    elif isinstance(mineralform, frozenset):
        bits = list(mineralform)
        random.shuffle(bits)
        for part in bits:
            quaffle = drawfrom(part, molmany)
            if quaffle:
                return quaffle
        return None
    else:
        bits = ()
        for part in mineralform:
            quaffle = drawfrom(part, molmany)
            if quaffle:
                bits = bits + quaffle
            else:
                for unroll in bits:
                    molmany[unroll] += 1
                return None
        return bits

def silicamelt(molmany):
    order = ['Na', 'K', 'Mg', 'Ca', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Al', 'Si', 'C', 'N', 'P', 'S', 'O', 'H', 'F', 'Cl', 'He', 'Ne', 'Ar']
    activeminerals = list(minerals)
    mineraldata = dict()

    while activeminerals:
        midx = random.randrange(0, len(activeminerals))
        mineral = activeminerals[midx]
        mineralname = mineral['name']
        
        components = drawfrom(mineral['form'], molmany)
        if components:
            if mineralname not in mineraldata:
                mineraldata[mineralname] = {'molecules': dict(),
                                            'forms':     dict()}
                
            atoms = dict()
            for mol in components:
                if mol in mineraldata[mineralname]['molecules']:
                    mineraldata[mineralname]['molecules'][mol] += 1
                else:
                    mineraldata[mineralname]['molecules'][mol] = 1
                    
                for atom, m in molecules[mol]['components'].items():
                    atoms[atom] = atoms[atom] + m if atom in atoms else m
                    
            formula = ''.join(atom + ('' if atoms[atom]==1 else str(atoms[atom])) for atom in order if atom in atoms)
            if formula in mineraldata[mineralname]['forms']:
                mineraldata[mineralname]['forms'][formula] += 1
            else:
                mineraldata[mineralname]['forms'][formula] = 1
        else:
            del activeminerals[midx]

    for m, d in mineraldata.items():
        d.update(process(''.join((f * q for f, q in d['forms'].items()))))
        scont = silica_content(d)
        sclass = silica_classify(scont)
        d['scont'] = scont
        d['sclass']   = sclass

    return mineraldata

if __name__ == '__main__':
    for mol in mollist:
        print(f'{mol} => {molecules[mol]!r}')


    molmany = {'FeO':10, 'MgO': 10, 'Na2O': 10, 'NaOH': 10, 'MgO2H2': 10, 'Al2O3':10, 'SiO2': 100}
    mineraldata=silicamelt(molmany)

    for m, d in mineraldata.items():
        print(f"{m:12s} {d['sclass']:>12s} ({d['scont']*100:6.2f}%)")
        print(f"  FORMS: {d['forms']!r}")
        print(f"  MOLES: {d['molecules']!r}")
        
    print(molmany)
    
