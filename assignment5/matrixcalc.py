import numpy as np
import itertools

s = np.array([[-1,0,0,0,0,0,0,1,0,0,0],
              [1,-1,0,-1,0,0,0,0,0,0,0],
              [0,1,-1,0,0,0,0,0,0,0,0],
              [0,0,1,1,-1,0,0,0,0,0,0],
              [0,0,0,0,1,-1,0,0,0,0,0],
              [0,0,0,0,0,1,-1,0,0,0,0],
              [0,0,0,0,0,0,1,0,0,0,-1],
              [-1,0,0,0,0,0,1,0,0,0,0],
              [1,0,0,0,0,0,-1,0,0,0,0],
              [0,0,0,-1,0,-1,0,0,2,0,0],
              [0,0,0,1,0,1,0,0,0,-2,0]])

if __name__ == '__main__':
    vecspace = [0,1,2]
    poss = itertools.product(vecspace, repeat=11)

    while True:
        try:
            nxt = np.array(next(poss))
            test = s @ nxt
            if not np.any(test):
                print(nxt)
        except:
            break

