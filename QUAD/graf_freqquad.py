import numpy as np 
import matplotlib.pyplot as plt 
from scipy.special import jn_zeros 

fig, ax = plt.subplots()
dades = np.loadtxt('QUAD/EX/quad_ex_freq.txt')
dadesT = np.transpose(dades)

k_llista = [] #K=W/PI TEORICS OBTINGUTS
deg = [] #DEGENERACIÓ DELS VALORS K

for m in range(0,40,2):
    for n in range(0,40,2):
        k = np.sqrt(n**2+m**2)
        if k not in k_llista:
            k_llista.append(k)
            deg.append(1)
        else:
            deg[k_llista.index(k)] += 1

sorted_k = np.sort(np.array(k_llista))

#CALCUL ERROR MIG DE LES FREQUENCIES RESONANTS/PI:#########################################################################################3
resonancia = [] #GUARDEM W,P DELS MAXIMS LOCALS

for i in range(1,len(dadesT[0])-1): #CONDICIÓ DE MAXIM
    if dadesT[2][i]> dadesT[2][i-1] and dadesT[2][i]> dadesT[2][i+1]:
        resonancia.append([dadesT[0][i], dadesT[2][i]])

resonanciaT = np.transpose(np.array(resonancia))

err = 0 #CONTADOR PER A L'ERROR
n = 0  
for res in resonanciaT[0]:
    for k in k_llista:
        if abs(k - res/np.pi)<0.15:
            err += abs(k - res/np.pi)
            n += 1
            
print(err/n)


#PLOT#########################################################################################################################################
ax.plot(dadesT[0]/np.pi, dadesT[2]/max(dadesT[2][:180]), label = r'Placa quadrada')

#representem els valors teorics k per a les maximes degeneracions (4,6)
k_llista = np.array(k_llista)
k_maxdeg = (k_llista[np.where(np.array(deg)==4)])
k_maxdeg = np.append(k_maxdeg, k_llista[np.where(np.array(deg)==6)])  

for i in range(1, len(resonanciaT[1])): #plot de valors teorics per k
        if sorted_k[i] <16.5:
            y = dadesT[2][np.where(np.abs(dadesT[0]/np.pi-sorted_k[i])<0.1)][-1]/np.max(dadesT[2][:180])
            plt.vlines(sorted_k[i],0,1.1, lw=1, alpha=0.75, color = 'lightblue')
  
for k in k_maxdeg:
    ax.vlines(k,0,1.25, lw = 2, color='lightblue', alpha =1)

ax.text(9.75, 1, r'$K_{0,10}$')
ax.text(15.9245155-0.05,1, r'$K_{2,16}$')
ax.set_xlim(0,17)
ax.set_ylim(0,1.1)
ax.set_xlabel(r'$\frac{\hat{k}}{\pi}$', fontsize=16)
ax.set_ylabel(r'$\frac{E}{E_{max}}$', fontsize=16)
plt.show()