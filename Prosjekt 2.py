import random as random
import numpy as np
import matplotlib.pyplot as plt


alfa = 1.0610
beta = 53.0516
a = alfa/12.
b = -alfa/24.
c = 0
d = -alfa/beta
N = 40

fjærkonst = 0.1
l = 5*10**(-9)
r = 10**(-9)
E = 300*10**6
tmax = 1

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


# Oppgave 2a:

def alpha(masse, l, E, r):
    return ((9.81*masse*l**3)/(E*3.14*r**4/4))
    
def beta_k(fjærkonst, l, E, r):
    return 10**(8)#((fjærkonst*l**3)/(E*3.14*r**4/4))

def u_vektor(N):
    u = []
    for i in range(N):
        u.append(-alfa/6.)
        u.append(-alfa/4.)
        u.append(-alfa/6.)
        u.append(-alfa/24.)
    return u

    
# def A_matrise(N, beta, alfa):
#     s = 4*N
#     A = np.zeros((s,s))
#     temp = 0
#     for i in range(0,s,4):
#         for j in range(i, i+4):
#             A[j][j] = 1
#         A[i][i+3] = beta[temp]/6
#         temp += 1
#         if i==0:
#             for j in range(4):
#                 A[i+j][s-(4-j)] = -1
#             A[i+1][s-4] = -3
#             A[i+2][s-4] = -3
#             A[i+2][s-3] = -2
#             for j in range(s-4,s):
#                 A[i+3][j] = -1
#         else:
#             for j in range(4):
#                 A[i+j][i-(4-j)] = -1
#             A[i+1][i-4] = -3
#             A[i+2][i-4] = -3
#             A[i+2][i-3] = -2
#             for j in range(i-4, i):
#                 A[i+3][j] = -1
#     return A

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

u = u_vektor(N)
beta_liste = np.zeros(N)
for i in range(N):
    beta_liste[i] = beta #beta_k(fjærkonst, l, E, r)
A = A_matrise(N, beta_liste, alfa)



ksi=np.linspace(0,N,1000*N)
eta_verdier=np.zeros(1000*N)
eta_derivert_verdier=np.zeros(1000*N)
eta_dobbeltderivert_verdier=np.zeros(1000*N)
eta_trippelderivert_verdier=np.zeros(1000*N)
konstanter = np.linalg.solve(A,u)

def eta(k,ksi, alfa, konstanter):
    return ((-alfa/24.)*((ksi-k)**4)+(konstanter[4*k]*(ksi-k)**3)
            +(konstanter[4*k+1]*(ksi-k)**2)+(konstanter[4*k+2]*(ksi-k))*konstanter[4*k+3])

def eta_derivert(k,ksi, alfa, konstanter):
    return ((-alfa/6.)*((ksi-k)**3)+(3*konstanter[4*k]*(ksi-k)**2)
            +(2*konstanter[4*k+1]*(ksi-k))+(konstanter[4*k+2]))

def eta_dobbeltderivert(k,ksi, alfa, konstanter):
    return ((-alfa/2.)*((ksi-k)**2)+(6*konstanter[4*k]*(ksi-k))
            +(2*konstanter[4*k+1]))

def eta_trippelderivert(k,ksi, alfa, konstanter):
    return (-(alfa*(ksi-k))+6*konstanter[4*k])

for i in range (N):
    for j in range (int(len(ksi)/N)):
        eta_verdier[j+i*len(ksi)/N]=eta(i,ksi[j+i*len(ksi)/N], alfa, konstanter)

for i in range(N):
    for j in range(int(len(ksi)/N)):
        eta_derivert_verdier[j+i*len(ksi)/N] = eta_derivert(i,ksi[j+i*len(ksi)/N], alfa, konstanter)

for i in range(N):
    for j in range(int(len(ksi)/N)):
        eta_dobbeltderivert_verdier[j+i*len(ksi)/N] = eta_dobbeltderivert(i,ksi[j+i*len(ksi)/N], alfa, konstanter)

for i in range(N):
    for j in range(int(len(ksi)/N)):
        eta_trippelderivert_verdier[j+i*len(ksi)/N] = eta_trippelderivert(i,ksi[j+i*len(ksi)/N], alfa, konstanter)





# plt.plot(ksi,eta_verdier)
# plt.savefig("Oppgave 2a,uderivert.pdf")
# plt.show()
# plt.plot(ksi,eta_derivert_verdier)
# plt.savefig("Oppgave 2a,derivert.pdf")
# plt.show()
# plt.plot(ksi,eta_dobbeltderivert_verdier)
# plt.savefig("Oppgave 2a,dobbeltderivert.pdf")
# plt.show()
# plt.plot(ksi,eta_trippelderivert_verdier,'bo')
# plt.savefig("Oppgave 2a,trippelderivert.pdf")
# plt.show()





def t_k(t_max, N):
    tk = np.zeros(N)
    for i in range(N):
        tk[i] = random.random() * t_max
        while tk[i] == 0:
            tk[i] = random.random() * t_max
    return tk
    
def r_k(beta, konstanter, tk, N):
    d_j = np.zeros(N)
    a = 0
    for i in range(len(konstanter)):
        if i%4 == 3:
            d_j[a]=(konstanter[i])
            a += 1
    rk = np.zeros(N)
    for i in range(N):
        rk[i] = (-beta[i] * d_j[i])/tk[i]
    return rk
    
    
def oppgave3(t_max, N):
    beta = np.zeros(N)
    for i in range(N):
        beta[i] = beta_k(fjærkonst, l, E, r)
    tk = t_k(t_max, N)
    #print(tk)
    
    ryker = np.zeros(N)
    alphaer = np.zeros(N)
    for i in range(N):
        A = A_matrise(N, beta, 1)
        u = u_vektor(N)
        konstanter = np.linalg.solve(A,u)
        rk = r_k(beta, konstanter, tk, N)
        temp = np.argmax(rk)
        alphaer[i] = beta[temp]/rk[temp]
        beta[temp] = 0
        ryker[i] = temp
    #print(ryker)
        
    return alphaer 

alphaer = oppgave3(tmax, N)
print(alphaer)

n_N_liste = np.zeros(N)
for i in range(N):
    n_N_liste[i]= i/N


def midlere_alpha(N):
    alpha_middel = oppgave3(tmax, N)
    for i in range(N):
        alpha2 = oppgave3(tmax, N)
        for j in range(N):
            alpha_middel[j] = (alpha_middel[j] + alpha2[j])
    return alpha_middel/N
    
plt.plot(n_N_liste, midlere_alpha(N))
plt.show()
            
            
        

