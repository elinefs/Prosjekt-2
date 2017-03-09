import numpy as np
import random as random
import Matrise as ma
import Konstanter as ko


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
    
    
def losrivelse(t_max, N):
    tk = t_k(t_max, N)
    #print(tk)
    beta = ko.beta_k(ko.fjaerkonst, ko.l, ko.E, ko.r)
    
    ryker = np.zeros(N)
    alphaer = np.zeros(N)
    for i in range(N):
        A = ma.A_matrise(N, beta)
        u = ma.u_vektor(N, 1)
        konstanter = np.linalg.solve(A,u)
        rk = r_k(beta, konstanter, tk, N)
        temp = np.argmax(rk)
        alphaer[i] = beta[temp]/rk[temp]
        beta[temp] = 0
        ryker[i] = temp
    #print(ryker)
        
    return alphaer 



def midlere_alpha(tmax, N):
    alpha_middel = losrivelse(tmax, N)
    for i in range(N):
        alpha2 = losrivelse(tmax, N)
        for j in range(N):
            alpha_middel[j] = (alpha_middel[j] + alpha2[j])
    return alpha_middel/N