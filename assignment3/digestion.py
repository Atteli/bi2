import argparse
import sys
import math

import enzymelib as ezlib
import fastareader as fr
import peptidemasshandler as ph

def digest(fasta, enz):
    enzyme_instance = ezlib.Enzyme(enz)
    return enzyme_instance.cleave(fasta)

def make_tuples(list1, list2, list3):
    print(list1, list2, list3)
    assert(len(list1) == len(list2) == len(list3))
    tpls = []
    for index, elem in enumerate(list1):
        tpls.append((list1[index], list2[index], list3[index]))
    return tpls

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="In Silico Digestor")
    parser.add_argument("--fasta", help="FASTA-Sequence")
    parser.add_argument("--enzyme", nargs='?', help="Enzyme for sequence digestion", type=str, default='trypsin')
    parser.add_argument("--masses", nargs='?', help="Masses of amino acids", type=str, default='./aa_masses_extended.csv')

    args = parser.parse_args()

    if args.fasta is None:
        print('Usage: python3 digestion.py --fasta {fastafile} [--enzyme {str}]')
        sys.exit('Error! No FASTA file specified. Exiting...')

    fasta = args.fasta
    enz = args.enzyme
    ms = args.masses

    headers, seqs = fr.fastaread(fasta)

    for i,s in enumerate(seqs):
        print("processing seq ", i+1 , " of ", len(seqs), ", of file ", fasta, "\n")
        scatter = digest(seqs[0], enz)
        mono_masses, avg_masses, total_mass_mono, total_mass_avg = ph.get_masses(ms, scatter, 'M')
        print("performing " + enz + " digestion of seq " + str(headers[i])[1:])
        print("total peptide mass (monoisotopic):\n", total_mass_mono)
        print("total peptide mass (average):\n", total_mass_avg)
        print("peptide pieces after digestion:\n", scatter)
        print("monoisotopic masses of pieces:\n", mono_masses)
        print("average masses of pieces:\n", avg_masses)


        print("sorted scatter list: (piece, monoisotopic_mass, avg_mass), sorted by piece length:")
        lst = make_tuples(scatter, mono_masses, avg_masses)
        print("", sorted(lst, key=lambda lst: len(lst[0])))


