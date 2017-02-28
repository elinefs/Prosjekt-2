__author__ = 'HaakonGryvill'

alfa = 1.
beta = 20.
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
"""u = np.transpose(u)
    for i in range(len(A)):
    print(len(A[i]),i)
    
    """

print(np.linalg.solve(A,u))

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
    for i in range(N):
        
        
        return A
