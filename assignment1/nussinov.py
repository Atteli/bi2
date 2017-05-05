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

def resultformat(mode='bpseq', verbose=False):
    if __SEQ and __PAIRS:
        temp = list(__SEQ)
        for pair in __PAIRS:
            if mode == 'bpseq':
                if verbose:
                    print(pair[0] + 1, __SEQ[pair[0]], __SEQ[pair[1]], pair[1] + 1)
            temp[pair[0]] = '('
            temp[pair[1]] = ')'
        if mode == 'dot-bracket':
            for i in range(len(temp)):
                if not (temp[i] == '(' or temp[i] == ')'):
                    temp[i] = '.'
            if verbose:
                print(''.join(temp))
            return ''.join(temp)
        return __PAIRS

def printresults(verbose=False):
    if verbose:
        print('\nSequence:', __SEQ)
        print('\n')

        print('bpseq:\n')
    bpseq = resultformat(mode='bpseq')

    if verbose:
        print('\ndot-bracket:\n')
    db = resultformat(mode='dot-bracket')

    if verbose:
        print('\n')
        print('Score:', __M[0, len(__M[0]) - 1])

    return bpseq, db

def onconsole():
    parser = argparse.ArgumentParser(description="Nussinov algorithm")
    parser.add_argument("--fasta", help="Query sequence fasta")
    parser.add_argument("--min_loop", nargs='?', help="Minimum loop length, default: 3", type=int, default=3)
    parser.add_argument("--GC", nargs='?', help="Score for GC / CG base pair, default: 0", type=int, default=0)
    parser.add_argument("--AU", nargs='?', help="Score for AU / UA base pair, default: 0", type=int, default=0)
    parser.add_argument("--GU", nargs='?', help="Score for GU / UG base pair, default: 0", type=int, default=0)

    args = parser.parse_args()
    print(args)
    if args.fasta is None:
        print(
            'Usage: python3 nussinov.py --fasta {fastafile} [--min_loop {int}] [--GC {int}] [--AU {int}] [--GU {int}]')
        sys.exit('Error! No FASTA file specified. Exiting...')

    headers, seqs = fastaread(args.fasta)

    global __SEQ
    global __GC
    global __AU
    global __GU
    global __MIN_LOOP

    __SEQ = seqs[0].upper().replace('T', 'U')
    __GC = args.GC
    __AU = args.AU
    __GU = args.GU
    __MIN_LOOP = args.min_loop

def onide(gc, au, gu, min_loop):
    headers, seqs = fastaread('E:/GitProjects/bi2/assignment1/supplement/at.fasta')

    global __SEQ
    global __GC
    global __AU
    global __GU
    global __MIN_LOOP

    __SEQ = seqs[0].upper().replace('T', 'U')
    __GC = gc
    __AU = au
    __GU = gu
    __MIN_LOOP = min_loop

def resetglobals():
    global __SEQ
    global __GC
    global __AU
    global __GU
    global __MIN_LOOP
    global __M
    global __PAIRS

    __GC = 0
    __AU = 0
    __GU = 0
    __MIN_LOOP = 0
    __M = []
    __SEQ = ''
    __PAIRS = []

def bruteforce(r):
    global __M
    cmpseq = '(((((((..((((.......)))).(((((.......))))).....((.((.......)).))))))))).'
    tmp = 0
    tmpdb = ''
    tmpscr = 0
    tmpargs = []
    for gc in range(r):
        for au in range(r):
            for gu in range(r):
                for min_loop in range(5,10):
                    resetglobals()
                    onide(gc, au, gu, min_loop)
                    __M = nussinovfill()
                    nussinovtraceback(0, len(__SEQ) - 1)
                    bpseq, db = printresults()
                    print(db, tmp / (r*r*r), __GC, __AU, __GU, __MIN_LOOP, __M[0, len(__M[0])-1])
                    if __M[0, len(__M[0])-1] > tmpscr:
                        tmpdb = db
                        tmpscr = __M[0, len(__M[0])-1]
                        tmpargs = [gc, au, gu, min_loop]
                    tmp += 1
                    if db == cmpseq:
                        print('found config!', __GC, __AU, __GU, __MIN_LOOP)
                        sys.exit()
    print('finished.')
    print(tmpdb, tmpscr, tmpargs)

if __name__ == '__main__':
    '''
    #onconsole()
    onide()

    __M = nussinovfill()
    nussinovtraceback(0, len(__SEQ) - 1)

    bpseq, db = printresults()
    '''
    bruteforce(10)




