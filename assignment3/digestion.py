import argparse
import sys
import pandas as pd

import assignment3.enzymelib as ezlib

def fastaread(fpath):
    '''
    :param fpath: File pach
    :return: list of headers / list of sequences of FASTA-file
    '''
    f = open(fpath, 'r')
    headers = []
    seqs = []
    temp = ''
    for line in f:
        if line:
            if line[0] == '>':
                headers.append(line)
                if temp:
                    seqs.append(temp)
                    temp = ''
            else:
                temp += line.strip()
    seqs.append(temp.upper())
    return headers, seqs

def get_masses(masses_sheet, matr):
    table = pd.read_csv(masses_sheet, sep=',')
    table = table.set_index("code")
    monoisotopic_masses = []
    average_masses = []
    for peptide_piece in matr:
        mono = 0
        avg = 0
        for c in peptide_piece:
            mono += table.loc[c].monoisotopic
            avg += table.loc[c].average
        monoisotopic_masses.append(mono)
        average_masses.append(avg)
    return monoisotopic_masses, average_masses

def digest(fasta, enz):
    enzyme_instance = ezlib.Enzyme(enz)
    return enzyme_instance.cleave(fasta)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="In Silico Digestor")
    parser.add_argument("--fasta", help="FASTA-Sequence")
    parser.add_argument("--enzyme", nargs='?', help="Enzyme for sequence digestion", type=str, default='trp')
    parser.add_argument("--masses", nargs='?', help="Masses of amino acids", type=str, default='./aa_masses.csv')

    args = parser.parse_args()

    if args.fasta is None:
        print('Usage: python3 digestion.py --fasta {fastafile} [--enzyme {str}]')
        sys.exit('Error! No FASTA file specified. Exiting...')

    fasta = args.fasta
    enz = args.enzyme
    ms = args.masses

    headers, seqs = fastaread(fasta)

    scatter = digest(seqs[0], enz)
    mono_masses, avg_masses = get_masses(ms, scatter)

    print(scatter)
    print(mono_masses)
    print(avg_masses)


