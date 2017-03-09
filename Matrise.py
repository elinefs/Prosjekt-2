import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve

def u_vektor(N, alfa):
    u = []
    for i in range(N):
        u.append(-alfa/6.)
        u.append(-alfa/4.)
        u.append(-alfa/6.)
        u.append(-alfa/24.)
    return u


def A_matrise(N, beta):
    s = 4*N
    A = np.identity(s)
    block = np.matrix('-1 0 0 0; -3 -1 0 0; -3 -2 -1 0;-1 -1 -1 -1')
    temp = 1
    for i in range(4,s,4):
        A[i:i+4,i-4:i] = block
        A[i][i+3] = beta[temp]/6
        temp += 1
    A[0:4,s-4:s] = block
    A[0][3] = beta[0]/6
    #A = sps.dok_matrix(A)
    return A
    
