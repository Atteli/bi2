from prody import *
import numpy as np
import os
import warnings


warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

confProDy(verbosity='none')

__RHALPHAHELICES = 0 #total alpha helox length
__RH310HELICES = 0 #total 3_10 helix length
__SHEETS = 0 #total beta sheets length
__TOTAL = 0 #total amino acids

__aminoacids = ['ALA', 'ARG',	'ASN', 'ASP', 'CYS', 'GLU',	'GLN', 'GLY', 'HIS', 'ILE',
                  'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

if __name__ == '__main__':
    '''
    parse all files in ./pdbs/ and extract lengths from alpha helices, 3_10 helices and beta sheets with ProDy package.
    Then, count all amino acids in all files.
    Finally, compute the amount of amino acids that take part in alpha helices, 3_10 helices and beta sheets.
    '''
    print('Depending on the amount of files to process, this may take a while...')
    for root, dirs, files in os.walk("./pdbs/", topdown=True):
        for name in files:
            atoms = parsePDB(os.path.join(root, name))
            header = parsePDB(os.path.join(root, name), model=0, header=True)
            for helix in header['helix_range']:
                if helix[2] == 1:
                    __RHALPHAHELICES += helix[5] - helix[4] + 1
                elif helix[2] == 5:
                    __RH310HELICES += helix[5] - helix[4] + 1
            for sheet in header['sheet_range']:
                __SHEETS += helix[5] - helix[4] + 1
            __TOTAL += len(list(filter(lambda item: item not in __aminoacids, atoms.getResnames())))
    #__TOTAL = __RHALPHAHELICES + __RH310HELICES + __SHEETS

    print('Alpha Helix Percentage:', (__RHALPHAHELICES / __TOTAL) * 100)
    print('3_10 Helix Percentage:',  (__RH310HELICES / __TOTAL) * 100)
    print('Beta Sheet Percentage:', (__SHEETS / __TOTAL) * 100)
