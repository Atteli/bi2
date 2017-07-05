import argparse
import sys
import numpy as np

import fastareader as fr
import peptidemasshandler as ph

def split_seq(seq):
    bmatr = []
    ymatr = []
    for i,e in enumerate(seq):
        bmatr.append(seq[:i+1])
        ymatr.append(seq[i:])
    return bmatr, ymatr

def calc_b_ion_masses(bmasses):
    return np.array(bmasses) + 1

def calc_y_ion_masses(ymasses):
    return np.array(ymasses)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MS experiment")
    parser.add_argument("--fasta", help="FASTA-Sequence")
    parser.add_argument("--masses", nargs='?', help="Masses of amino acids", type=str, default='./aa_masses.csv')

    args = parser.parse_args()

    if args.fasta is None:
        print('Usage: python3 digestion.py --fasta {fastafile} [--enzyme {str}]')
        sys.exit('Error! No FASTA file specified. Exiting...')

    fasta = args.fasta
    ms = args.masses

    headers, seqs = fr.fastaread(fasta)

    bmatr, ymatr = split_seq(seqs[0])

    mode = 'MH+'
    b_mono_masses, b_avg_masses, b_total_mass_mono, b_total_mass_avg = ph.get_masses(ms, bmatr)
    y_mono_masses, y_avg_masses, y_total_mass_mono, y_total_mass_avg = ph.get_masses(ms, ymatr, mode)

    print("b ions: ", bmatr)
    print("y ions:", ymatr)
    print("b ion masses: ", calc_b_ion_masses(b_mono_masses))
    print("y ion masses: ", calc_y_ion_masses(y_mono_masses))

