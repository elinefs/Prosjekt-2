
alfa = 1.0610
beta = 53.0516
a = alfa/12.
b = -alfa/24.
c = 0
d = -alfa/beta
N = 4
import numpy as np
import matplotlib.pyplot as plt
#hudw
# Oppgave 1:
A1 = [[1,0,0,beta/6,0,0,0,0,0,0,0,0,-1,0,0,0], # k = 0
     [0,1,0,0,0,0,0,0,0,0,0,0,-3,-1,0,0],
     [0,0,1,0,0,0,0,0,0,0,0,0,-3,-2,-1,0],
     [0,0,0,1,0,0,0,0,0,0,0,0,-1,-1,-1,-1],
     [-1,0,0,0,1,0,0,beta/6,0,0,0,0,0,0,0,0],# k = 1
     [-3,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
     [-3,-2,-1,0,0,0,1,0,0,0,0,0,0,0,0,0],
     [-1,-1,-1,-1,0,0,0,1,0,0,0,0,0,0,0,0],
     [0,0,0,0,-1,0,0,0,1,0,0,beta/6,0,0,0,0], # k = 2
     [0,0,0,0,-3,-1,0,0,0,1,0,0,0,0,0,0],
     [0,0,0,0,-3,-2,-1,0,0,0,1,0,0,0,0,0],
     [0,0,0,0,-1,-1,-1,-1,0,0,0,1,0,0,0,0],
     [0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0,beta/6], # k = 3
     [0,0,0,0,0,0,0,0,-3,-1,0,0,0,1,0,0],
     [0,0,0,0,0,0,0,0,-3,-2,-1,0,0,0,1,0],
     [0,0,0,0,0,0,0,0,-1,-1,-1,-1,0,0,0,1]]

u1 = [-alfa/6,-alfa/4.,-alfa/6.,-alfa/24.,-alfa/6,-alfa/4,-alfa/6,-alfa/24,
     -alfa/6,-alfa/4,-alfa/6,-alfa/24,-alfa/6,-alfa/4,-alfa/6,-alfa/24]


print(np.linalg.solve(A1,u1))

# Oppgave 2a:

def u_vektor(N):
    u = []
    for i in range(N):
        u.append(-alfa/6.)
        u.append(-alfa/4.)
        u.append(-alfa/6.)
        u.append(-alfa/24.)
    return u

    
def A_matrise(N):
    s = 4*N
    A = np.zeros((s,s))
    for i in range(0,s,4):
        for j in range(i, i+4):
            A[j][j] = 1
        A[i][i+3] = beta/6#(i/4)
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

u = u_vektor(N)
A = A_matrise(N)




def beta(i):
    return float(kappa(k)*(l**3)/B)


ksi=np.linspace(0,N,1000*N)
eta_verdier=np.zeros(1000*N)
eta_derivert_verdier=np.zeros(1000*N)
eta_dobbeltderivert_verdier=np.zeros(1000*N)
eta_trippelderivert_verdier=np.zeros(1000*N)
konstanter = np.linalg.solve(A,u)

def eta(k,ksi):
    return ((-alfa/24.)*((ksi-k)**4)+(konstanter[4*k]*(ksi-k)**3)
            +(konstanter[4*k+1]*(ksi-k)**2)+(konstanter[4*k+2]*(ksi-k))*konstanter[4*k+3])

def eta_derivert(k,ksi):
    return ((-alfa/6.)*((ksi-k)**3)+(3*konstanter[4*k]*(ksi-k)**2)
            +(2*konstanter[4*k+1]*(ksi-k))+(konstanter[4*k+2]))

def eta_dobbeltderivert(k,ksi):
    return ((-alfa/2.)*((ksi-k)**2)+(6*konstanter[4*k]*(ksi-k))
            +(2*konstanter[4*k+1]))

def eta_trippelderivert(k,ksi):
    return (-(alfa*(ksi-k))+6*konstanter[4*k])

for i in range (N):
    for j in range (len(ksi)/N):
        eta_verdier[j+i*len(ksi)/N]=eta(i,ksi[j+i*len(ksi)/N])

for i in range(N):
    for j in range(len(ksi)/N):
        eta_derivert_verdier[j+i*len(ksi)/N] = eta_derivert(i,ksi[j+i*len(ksi)/N])

for i in range(N):
    for j in range(len(ksi)/N):
        eta_dobbeltderivert_verdier[j+i*len(ksi)/N] = eta_dobbeltderivert(i,ksi[j+i*len(ksi)/N])

for i in range(N):
    for j in range(len(ksi)/N):
        eta_trippelderivert_verdier[j+i*len(ksi)/N] = eta_trippelderivert(i,ksi[j+i*len(ksi)/N])





plt.plot(ksi,eta_verdier)
plt.savefig("Oppgave 2a,uderivert.pdf")
plt.show()
plt.plot(ksi,eta_derivert_verdier)
plt.savefig("Oppgave 2a,derivert.pdf")
plt.show()
plt.plot(ksi,eta_dobbeltderivert_verdier)
plt.savefig("Oppgave 2a,dobbeltderivert.pdf")
plt.show()
plt.plot(ksi,eta_trippelderivert_verdier,'bo')
plt.savefig("Oppgave 2a,trippelderivert.pdf")
plt.show()
