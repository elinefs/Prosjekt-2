import numpy as np
import matplotlib.pyplot as plt

import Matrise as ma
import Konstanter as ko
import Løsrivelse as lo



def eta(k,ksi, alfa, konstanter):
    return ((-alfa/24.)*((ksi-k)**4)+(konstanter[4*k]*(ksi-k)**3)
            +(konstanter[4*k+1]*(ksi-k)**2)+(konstanter[4*k+2]*(ksi-k))+konstanter[4*k+3])

def eta_analytisk(k,ksi,alfa, beta):
    return ((-alfa/24.)*((ksi-k)**4)+(ko.a*(ksi-k)**3)
            +(ko.b*(ksi-k)**2)+(ko.c*(ksi-k))+ko.d)

def eta_derivert(k,ksi, alfa, konstanter):
    return ((-alfa/6.)*((ksi-k)**3)+(3*konstanter[4*k]*(ksi-k)**2)
            +(2*konstanter[4*k+1]*(ksi-k))+(konstanter[4*k+2]))

def eta_dobbeltderivert(k,ksi, alfa, konstanter):
    return ((-alfa/2.)*((ksi-k)**2)+(6*konstanter[4*k]*(ksi-k))
            +(2*konstanter[4*k+1]))

def eta_trippelderivert(k,ksi, alfa, konstanter):
    return (-(alfa*(ksi-k))+6*konstanter[4*k])
    
    
def plott1():
    N = ko.N
    beta = ko.beta
    alfa = ko.alfa
    A = ma.A_matrise(N, beta)
    u = ma.u_vektor(N, alfa)
    
    ksi = np.linspace(0,N,1000*N)
    eta_verdier=np.zeros(1000*N)
    eta_verdier_analytisk= np.zeros(1000*N)
    eta_derivert_verdier=np.zeros(1000*N)
    eta_dobbeltderivert_verdier=np.zeros(1000*N)
    eta_trippelderivert_verdier=np.zeros(1000*N)
    konstanter = np.linalg.solve(A,u)

    
    for i in range (N):
        for j in range (int(len(ksi)/N)):
            eta_verdier[j+i*len(ksi)/N]=eta(i,ksi[j+i*len(ksi)/N], alfa, konstanter)
    
    for i in range(N):
        for j in range(int(len(ksi)/N)):
            eta_verdier_analytisk[j+i*len(ksi)/N] = eta_analytisk(i,ksi[j+i*len(ksi)/N],alfa,beta)
    
    for i in range(N):
        for j in range(int(len(ksi)/N)):
            eta_derivert_verdier[j+i*len(ksi)/N] = eta_derivert(i,ksi[j+i*len(ksi)/N], alfa, konstanter)
    
    for i in range(N):
        for j in range(int(len(ksi)/N)):
            eta_dobbeltderivert_verdier[j+i*len(ksi)/N] = eta_dobbeltderivert(i,ksi[j+i*len(ksi)/N], alfa, konstanter)
    
    
    for i in range(N):
        for j in range(int(len(ksi)/N)):
            eta_trippelderivert_verdier[j+i*len(ksi)/N] = eta_trippelderivert(i,ksi[j+i*len(ksi)/N], alfa, konstanter)
    
    differanse = np.zeros(1000*N)
    for i in range(1000*N):
        differanse[i] = abs(eta_verdier[i]-eta_verdier_analytisk[i])
    
    plt.plot(ksi,differanse)
    plt.savefig("differanse.pdf")
    plt.show()
    
    
    
    plt.subplot(4,1,1)
    plt.plot(ksi,eta_verdier)
    
    plt.subplot(4,1,2)
    plt.plot(ksi,eta_derivert_verdier,'g')
    
    
    plt.subplot(4,1,3)
    plt.plot(ksi,eta_dobbeltderivert_verdier,'r')
    
    
    for i in range(N):
        plt.subplot(4,1,4)
        plt.plot(ksi[1000*i:1000*i+1000],eta_trippelderivert_verdier[1000*i:1000*i+1000],'k')
    
    
    plt.savefig("Oppgave2a.pdf")
    plt.show()
    
    plt.subplot(2,1,1)
    plt.plot(ksi,eta_verdier)
    
    plt.subplot(2,1,2)
    plt.plot(ksi,eta_verdier_analytisk,'g')
    
    plt.savefig("numerisk.pdf")
    plt.show()






def plott2():
    N = ko.N
    n_N_liste = np.zeros(N)
    for i in range(N):
        n_N_liste[i]= i/N
    
    midlere_alpha = lo.midlere_alpha(ko.tmax, ko.N)
        
    plt.plot(n_N_liste, midlere_alpha)
    plt.ylabel(r"$\alpha_n$")
    plt.xlabel("$n/N$")
    plt.show()