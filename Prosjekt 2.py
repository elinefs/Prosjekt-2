
alfa = 1.
beta = 2.
a = alfa/12.
b = -alfa/24.
c = 0
d = -alfa/beta
import numpy as np

# Oppgave 1:
A = [[6,0,0,beta,0,0,0,0,0,0,0,0,-6,0,0,0], # k = 0
     [0,1,0,0,0,0,0,0,0,0,0,0,-3,-1,0,0],
     [0,0,1,0,0,0,0,0,0,0,0,0,-3,-2,-1,0],
     [0,0,0,1,0,0,0,0,0,0,0,0,-1,-1,-1,-1],
     [-6,0,0,0,6,0,0,beta,0,0,0,0,0,0,0,0],# k = 1
     [-3,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
     [-3,-2,-1,0,0,0,1,0,0,0,0,0,0,0,0,0],
     [-1,-1,-1,-1,0,0,0,1,0,0,0,0,0,0,0,0],
     [0,0,0,0,-6,0,0,0,6,0,0,beta,0,0,0,0], # k = 2
     [0,0,0,0,-3,-1,0,0,0,1,0,0,0,0,0,0],
     [0,0,0,0,-3,-2,-1,0,0,0,1,0,0,0,0,0],
     [0,0,0,0,-1,-1,-1,-1,0,0,0,1,0,0,0,0],
     [0,0,0,0,0,0,0,0,-6,0,0,0,6,0,0,beta], # k = 3
     [0,0,0,0,0,0,0,0,-3,-1,0,0,0,1,0,0],
     [0,0,0,0,0,0,0,0,-3,-2,-1,0,0,0,1,0],
     [0,0,0,0,0,0,0,0,-1,-1,-1,-1,0,0,0,1]]

u = [-alfa,-alfa/4.,-alfa/6.,-alfa/24.,-alfa,-alfa/4,-alfa/6,-alfa/24,
     -alfa,-alfa/4,-alfa/6,-alfa/24,-alfa,-alfa/4,-alfa/6,-alfa/24]


#print(np.linalg.solve(A,u))

# Oppgave 2a:
N = 10
u_test = [-alfa,-alfa/4.,-alfa/6.,-alfa/24.]
u = []
def u_vektor(N):
    for i in range(N):
        u.append(-alfa)
        u.append(-alfa/4.)
        u.append(alfa/6.)
        u.append(-alfa/24.)
    return u
    
    
def A_matrise(N):
    s = 4*N
    a = np.zeros((s,s))
    for i in range(0,s,4):
        for j in range(i, i+4):
            A[j][j] = 1
        A[i][i+3] = beta#(i/4)/6
        if i==0:
            for j in range(4):
                A[i+j][s-(4-j)] = -1
            A[i+1][s-4] = -3
            A[i+2][s-4] = -3
            A[i+2][s-3] = -2
            for j in range(s-4,s):
                A[i+3][j] = -1
        else:
            for j in range(4):
                A[i+j][i-(4-j)] = -1
            A[i+1][i-4] = -3
            A[i+2][i-4] = -3
            A[i+2][i-3] = -2
            for j in range(i-4, i):
                A[i+3][j] = -1
            
    return A


for i in range(16):
    A = A_matrise(4)
    print(A[i])