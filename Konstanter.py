import numpy as np

alfa = 1.0610
beta = 53.0516

a = alfa/12.
b = -alfa/24.
c = 0
d = -alfa/beta

N = 20

fjaerkonst = 0.1
l = 5*10**(-9)
r = 10**(-9)
E = 300*10**6
tmax = 1

    
def beta_k(fjaerkonst, l, E, r):
    beta_liste = np.zeros(N)
    for i in range(N):
        beta_liste[i] = ((fjaerkonst*l**3)/(E*np.pi*r**4/4))
    return beta_liste

beta = beta_k(fjaerkonst, l, E, r)