import numpy as np


__GC = 0
__AU = 0
__GU = 0
__MIN_LOOP = 0
__M = []
__SEQ = ''

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
            print (i + 1, __SEQ[i], __SEQ[j], j + 1)
            nussinovtraceback(i + 1, j - 1)
        else:
            for k in range (i + 1, j - 1):
                if __M[i, j] == __M[i, k] + __M[k + 1, j]:
                    nussinovtraceback(i, k)
                    nussinovtraceback(k + 1, j)
                    break

    return 0

if __name__ == '__main__':
    seq1 = 'GUUCAUAAGAGGUCAACAGCAACGGGUGUU'
    seq2 = 'GGGAAACCU'

    __SEQ = seq1
    __GC = 1
    __AU = 10
    __GU = 0
    __MIN_LOOP = 5

    __M = nussinovfill()
    print(__M)
    nussinovtraceback(0, len(__SEQ) - 1)
    print('Score:', __M[0, len(__M[0]) - 1])
