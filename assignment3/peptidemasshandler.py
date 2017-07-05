import pandas as pd
import math

def __truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n

def charge(table, mode):
    ch_m = 0
    ch_a = 0
    if mode == 'M':
        ch_m = (2 * table.loc['HYP'].monoisotopic) + table.loc['OXY'].monoisotopic
        ch_a = (2 * table.loc['HYP'].average) + table.loc['OXY'].average
    elif mode == 'MH+':
        ch_m = (3 * table.loc['HYP'].monoisotopic) + table.loc['OXY'].monoisotopic
        ch_a = (3 * table.loc['HYP'].average) + table.loc['OXY'].average
    return ch_m, ch_a

def get_masses(masses_sheet, matr, mode='null'):
    table = pd.read_csv(masses_sheet, sep=',')
    table = table.set_index("code")

    ch_m, ch_a = charge(table, mode)

    monoisotopic_masses = []
    average_masses = []
    fullseq = str(''.join(matr))
    tempmatr = list(matr)
    tempmatr.append(fullseq)
    for peptide_piece in tempmatr:
        mono = 0
        avg = 0
        for c in peptide_piece:
            mono += table.loc[c].monoisotopic
            avg += table.loc[c].average
        monoisotopic_masses.append(__truncate(mono + ch_m, 3))
        average_masses.append(__truncate(avg + ch_a, 3))

    return monoisotopic_masses[:-1], average_masses[:-1], \
           monoisotopic_masses[len(monoisotopic_masses)-1], average_masses[len(average_masses)-1]
