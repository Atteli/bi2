import numpy as np
import argparse
import sys


__GC = 0
__AU = 0
__GU = 0
__MIN_LOOP = 0
__M = []
__SEQ = ''
__PAIRS = []

def fastaread(fpath):
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
                temp += line
    seqs.append(temp)
    return headers, seqs

def deltafunc(i, j):
    test = ord(__SEQ[i]) + ord(__SEQ[j])
    return __GC if test == 138 else __AU if test == 150 else __GU if test == 156 else 0

def nussinovfill():
    m = np.zeros((len(__SEQ), len(__SEQ)))
    for l in range(1, len(__SEQ)):
        for j in range(l, len(__SEQ)):
            i = j - l
            m[i, j] = np.max([m[i + 1, j],
                              m[i, j - 1],
                              m[i + 1, j - 1] + deltafunc(i, j) if __MIN_LOOP < j - i else 0,
                              np.max(0 if j - i < 2 else [m[i, k] + m[k + 1, j] for k in range(i + 1, j)])])
    return m

def nussinovtraceback(i, j):
    if i < j:
        if __M[i, j] == __M[i + 1, j]:
            nussinovtraceback(i + 1, j)
        elif __M[i,j] == __M[i, j-1]:
            nussinovtraceback(i, j - 1)
        elif __M[i,j] == __M[i + 1, j - 1] + deltafunc(i, j):
            #print (i + 1, __SEQ[i], __SEQ[j], j + 1)
            __PAIRS.append([i, j])
            nussinovtraceback(i + 1, j - 1)
        else:
            for k in range (i + 1, j - 1):
                if __M[i, j] == __M[i, k] + __M[k + 1, j]:
                    nussinovtraceback(i, k)
                    nussinovtraceback(k + 1, j)
                    break

    return 0

def resultformat(mode='bpseq'):
    if __SEQ and __PAIRS:
        temp = list(__SEQ)
        for pair in __PAIRS:
            if mode == 'bpseq':
                print(pair[0] + 1, __SEQ[pair[0]], __SEQ[pair[1]], pair[1] + 1)
            temp[pair[0]] = '('
            temp[pair[1]] = ')'
        if mode == 'dot-bracket':
            for i in range(len(temp)):
                if not (temp[i] == '(' or temp[i] == ')'):
                    temp[i] = '.'
            print(''.join(temp))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Nussinov algorithm")
    parser.add_argument("--fasta", help="Query sequence fasta")
    parser.add_argument("--min_loop", nargs='?', help="Minimum loop length, default: 3", type=int, default=3)
    parser.add_argument("--GC", nargs='?', help="Score for GC / CG base pair, default: 0", type=int, default=0)
    parser.add_argument("--AU", nargs='?', help="Score for AU / UA base pair, default: 0", type=int, default=0)
    parser.add_argument("--GU", nargs='?', help="Score for GU / UG base pair, default: 0", type=int, default=0)

    args = parser.parse_args()
    print(args)
    if args.fasta is None:
        print('Usage: python3 nussinov.py --fasta {fastafile} [--min_loop {int}] [--GC {int}] [--AU {int}] [--GU {int}]')
        sys.exit('Error! No FASTA file specified. Exiting...')


    headers, seqs = fastaread(args.fasta)

    #print(headers, seqs)

    __SEQ = seqs[0]
    __GC = args.GC
    __AU = args.AU
    __GU = args.GU
    __MIN_LOOP = args.min_loop

    print('\nSequence:', __SEQ)
    print('\n')

    __M = nussinovfill()
    nussinovtraceback(0, len(__SEQ) - 1)

    print('bpseq:\n')
    resultformat(mode='bpseq')

    print('\ndot-bracket:\n')
    resultformat(mode='dot-bracket')

    print('\n')
    print('Score:', __M[0, len(__M[0]) - 1])

