from prody import *
import os

__SHEETS = 0
__RHALPHAHELICES = 0
__RH310HELICES = 0
__TOTAL = 0

if __name__ == '__main__':
    aminoacids = ['ALA', 'ARG',	'ASN', 'ASP', 'CYS', 'GLU',	'GLN', 'GLY', 'HIS', 'ILE',
                  'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    atoms = parsePDB('1a49')
    header = parsePDB('1a49', model=0, header=True)

    print(atoms)
    print(len(list(filter(lambda item: item not in aminoacids, atoms.getResnames()))))
    print(header['sheet_range'])


