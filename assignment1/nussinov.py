import numpy as np


def deltafunc(seq, i, j, gc, au, gu):
    test = ord(seq[i]) + ord(seq[j])
    return au if test == 150 else gc if test == 138 else gu if test == 156 else 0

def nussinovfill(seq:str, gc=2, au=3, gu=1, min_loop=0):
    m = np.zeros((len(seq), len(seq)))
    for l in range(1, len(seq)):
        for j in range(l, len(seq)):
            i = j - l
            m[i, j] = np.max([m[i + 1, j],
                              m[i, j - 1],
                              m[i + 1, j - 1] + deltafunc(seq, i, j, gc, au, gu) if min_loop < j - i else 0,
                              np.max(0 if j - i < 2 else [m[i, k] + m[k + 1, j] for k in range(i + 1, j)])])
    return m

def nussinovtraceback(seq, m, i, j, gc, au, gu):
    if i < j:
        if m[i, j] == m[i + 1, j]:
            nussinovtraceback(seq, m, i + 1, j, gc, au, gu)
        elif m[i,j] == m[i, j-1]:
            nussinovtraceback(seq, m, i, j - 1, gc, au, gu)
        elif m[i,j] == m[i + 1, j - 1] + deltafunc(seq, i, j, gc, au, gu):
            print (i + 1, seq[i], seq[j], j + 1)
            nussinovtraceback(seq, m, i + 1, j - 1, gc, au, gu)
        else:
            for k in range (i + 1, j - 1):
                if m[i, j] == m[i, k] + m[k + 1, j]:
                    nussinovtraceback(seq, m, i, k, gc, au, gu)
                    nussinovtraceback(seq, m, k + 1, j, gc, au, gu)
                    break

    return 0

if __name__ == '__main__':
    seq1 = 'GUUCAUAAGAGGUCAACAGCAACGGGUGUU'
    seq2 = 'GGGAAACCU'

    seq = seq1
    gc = 1
    au = 10
    gu = 0
    min_loop = 5

    m = nussinovfill(seq, gc, au, gu, min_loop)
    print(m)
    nussinovtraceback(seq, m, 0, len(seq) - 1, gc, au, gu)
    print('Score:', m[0, len(m[0]) - 1])
