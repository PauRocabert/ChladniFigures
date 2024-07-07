import numpy as np 
import matplotlib.pyplot as plt 
from scipy.special import jn_zeros 
from scipy.special import jv
from scipy.stats import linregress

plt.style.use('seaborn-darkgrid')

dades = np.loadtxt('MNII_PF/RAD/rad_freq.txt')
dadesT = np.transpose(dades)

a = jn_zeros(1,10) #W TEORICS

resonancia = [] #GUARDEM W,P DELS MAXIMS LOCALS
for i in range(1,len(dadesT[0])-1):
    if dadesT[2][i]> dadesT[2][i-1] and dadesT[2][i]> dadesT[2][i+1]:
        resonancia.append([dadesT[0][i], dadesT[2][i]])

resonanciaT = np.transpose(np.array(resonancia))
print(resonanciaT[0])

#CALCUL ERROR MIG DE LES FREQUENCIES RESONANTS/PI:#########################################################################################3
error = 0
n = 0
for alpha in a:
    for res in resonanciaT[0]:
        if abs(res - alpha) < 0.15:
            n+= 1
            error+= abs(res - alpha)
print(error/n)


#PLOT#######################################################################################
fig, ax = plt.subplots(figsize=(18,12))

for i in range(len(a)):
    ax.vlines(a[i],0,1.1, lw =1)
    ax.text(a[i]-0.1, 1.05, r'$\alpha^{1}'+f'_{i}$')

ax.plot(dadesT[0][10:] ,dadesT[2][10:]/np.max(dadesT[2][10:]), label='Placa circular')
ax.set_xlim(2,20)
ax.set_ylim(0,1.1)
ax.set_xlabel(r'$\frac{\hat{k}}{\pi}$', fontsize=22)
ax.set_ylabel(r'$\frac{E}{E_{max}}$', fontsize=22)

plt.show()
