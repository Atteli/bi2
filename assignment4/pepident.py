import argparse
import sys
import pandas as pd

import assignment3.peptidemasshandler as ph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MS experiment")
    parser.add_argument("--msspectrum", help="Tandem MS spectrum of an unknown peptide")
    parser.add_argument("--masses", nargs='?', help="Masses of amino acids", type=str, default='../assignment3/aa_masses_extended.csv')

    args = parser.parse_args()

    if args.msspectrum is None:
        print('Usage: python3 pepident.py --msspectrum {tandem ms spectrum file}')
        sys.exit('Error! No spectrum file specified. Exiting...')

    spectrum = args.msspectrum
    ms = args.masses

    table = pd.read_csv(ms, sep=',')
    table = table.set_index("monoisotopic")

    ch_m = (3 * 1.00728) + 15.9949146

    print(spectrum, table, ch_m)
