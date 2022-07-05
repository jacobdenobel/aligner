import numpy as np
from random import choice
from time import time
from gc import collect

class Align():
    
    cost = {"match":1,  "mismatch":-1, "gap":-2}
    aplhabet = "ATCG"
   
    def __init__(self, n, m):
        s1 = self.generateDNA(n)
        s2 = self.generateDNA(m)
        col = n + 1
        row = m + 1

        matrix = np.zeros((col,row), dtype=int) 
        matrix[:,0] = np.arange(col) * self.cost['gap']       
        matrix[0,:] = np.arange(row) * self.cost['gap']

        dmatrix = np.zeros((col,row), dtype=int) 
        dmatrix[:,0] = 1
        dmatrix[0,:] = 2
        dmatrix[0][0] = -1

        for i in range(1, col):
            for j in range(1, row):
                directions = (
                    matrix[i-1][j-1] + (1 if s1[i-1]==s2[j-1] else -1),
                    matrix[i-1][j] + -2,
                    matrix[i][j-1] + -2,
                )
                M = max(enumerate(directions), key=lambda l:l[1])[0]

                matrix[i][j] = directions[M]
                dmatrix[i][j]= M

        i=j=1
        o1, o2 ="", "",

        while dmatrix[col-i][row-j] != -1:
            current = dmatrix[col-i][row-j]
            if current == 0:
                i+=1
                j+=1
                o1 += s1[col-i]
                o2 += s2[row-j]
            elif current == 1:
                i+=1
                o1 += "-"
                o2 += s2[row-j]
            else:
                j+=1
                o1 += s1[col-i]
                o2 += "-"

            

        print("Score", matrix[col-1][row-1])
      
        # print(o1)
        # print(o2)

    def generateDNA(self, l):
        return ''.join(choice(self.aplhabet) for _ in range(l))
    


if __name__ == "__main__":
   
    
    for n in (10, 100, 1000, 10000):
        tic = time()
        m = n 
        Align(n, m)
        collect()
        print(n,"Time Elapsed", time() - tic, "\n")
