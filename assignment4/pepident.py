import argparse
import sys
import pandas as pd
import numpy as np

#identification of sequence
def ident_seq(spectrum, aa_table, ch_m_y, ch_m_b, acc):
    actual_b_mass = 0
    actual_y_mass = 0
    seq_b = []
    seq_y = []
    aa_table = aa_table.sort_values(by=['monoisotopic'], ascending=[True])
    aa_table.reset_index(drop=True, inplace=True)

    #go over spectrum values
    for ind1, mass in enumerate(spectrum):
        b_candidates = []
        y_candidates = []
        actual_b_am = 0
        actual_y_am = 0
        actual_b_offset = 1
        actual_y_offset = 1

        #test if an amino acid fits to spectrum value
        for ind2, am in enumerate(aa_table.monoisotopic):
            b_offset = np.abs((mass - ch_m_b - actual_b_mass - am))
            y_offset = np.abs((mass - ch_m_y - actual_y_mass - am))

            if b_offset < acc:
                b_candidates.append(aa_table.code[ind2])
                actual_b_am = am
                actual_b_offset = b_offset
            elif y_offset < acc:
                y_candidates.append(aa_table.code[ind2])
                actual_y_am = am
                actual_y_offset = y_offset
            else:
                continue

        if len(b_candidates) and actual_b_offset < actual_y_offset:
            actual_b_mass += actual_b_am
            seq_b.append(b_candidates)
        if len(y_candidates) and actual_y_offset < actual_b_offset:
            actual_y_mass += actual_y_am
            seq_y.append(y_candidates)

    return seq_b, seq_y

#tool to print out the sequence as string instead of list
def pretty_print(lst):
    result = []
    for item in lst:
        if len(item) > 1:
            result.extend('[')
            for ind, elem in enumerate(item):
                result.extend(elem)
                if ind < len(item) - 1:
                    result.extend(',')
            result.extend(']')
        else:
            result.extend(item)

    return ''.join(result)

if __name__ == "__main__":
    #arg parsing
    parser = argparse.ArgumentParser(description="MS experiment")
    parser.add_argument("--msspectrum", help="Tandem MS spectrum of an unknown peptide")
    parser.add_argument("--masses", nargs='?', help="Masses of amino acids", type=str, default='./aa_masses.csv')

    args = parser.parse_args()

    if args.msspectrum is None:
        print('Usage: python3 pepident.py --msspectrum {tandem ms spectrum file}')
        sys.exit('Error! No spectrum file specified. Exiting...')

    #file parsing and data preparation
    spectrumfile = list(open(args.msspectrum))
    aa_table = pd.read_csv(args.masses, sep=',')
    spectrum = np.array([spectrumfile[i].strip('\n') for i in range(len(spectrumfile))], dtype=np.float32)
    spectrum.sort()

    #mh+ charge
    ch_m_b = 1.00728
    ch_m_y = (3 * 1.00728) + 15.9949146

    #mass accuracy (Da); the absolute of +- 0.05 Da is 0.1 Da.
    acc = 0.1

    #result output
    bseq, yseq = ident_seq(spectrum, aa_table, ch_m_y, ch_m_b, acc)
    print('sequence (b ion construction):', pretty_print(bseq))
    print('sequence (y ion construction):', pretty_print(yseq))


