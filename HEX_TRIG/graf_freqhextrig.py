import numpy as np 
import matplotlib.pyplot as plt 

plt.style.use('seaborn-darkgrid')
fig, ax = plt.subplots(figsize=(15,10))

dades = np.loadtxt('MNII_PF/HEX_TRIG/HEX/hex_freq.txt')
dadesT = np.transpose(dades)

dades2 = np.loadtxt('MNII_PF/HEX_TRIG/TRIG/trig_freq.txt')
dadesT2 = np.transpose(dades2)

#Valors teorics per K/PI#####################################################################################################################
#HEXAGON/TRIANGLE##################     
k_llista = []
parells =[]
deg = []
for m in range(0,16,1):
    for n in range(0,16,1):
            k = 4/(3)*np.sqrt(n**2+(m**2)- n*m)
            if k not in k_llista:
                k_llista.append(k)
                parells.append([m,n])
                deg.append(1)
            else:
                deg[k_llista.index(k)] += 1
                
#CALCUL ERROR MIG DE LES FREQUENCIES RESONANTS/PI:#########################################################################################3
#HEXAGON##################     
resonancia = []

for i in range(1,len(dadesT[0])-1):
    if dadesT[2][i]> dadesT[2][i-1] and dadesT[2][i]> dadesT[2][i+1]:
        resonancia.append([dadesT[0][i], dadesT[2][i]])

resonanciaT = np.transpose(np.array(resonancia))       

error = 0 #contadors
n = 0
for res in resonanciaT[0]:
    for num in k_llista:
        if num<5.5: #calcular pels valors mostrats
            if abs(num - res/np.pi)<0.15:
                n += 1
                error+= abs(num - res/np.pi)

print("hexagon:"+str(error/n))

#TRIANGLE##################
error = 0 #contadors
n = 0
for res in resonanciaT[0]:
    for m in range(0,16,1):
        for n in range(0,16,1):
            if (m+n)%3==0:
                num = 4/(3)*np.sqrt(n**2+(m**2)- n*m)
                if num<6.5 : #calcular pels valors mostrats
                    if abs(num - res/np.pi)<0.35:
                        n += 1
                        error+= abs(num - res/np.pi)

print("triangle:"+str(error/n))

#PLOT#########################################################################################################################################

for i in range(len(k_llista)):
    if k_llista[i] > 1.0 and k_llista[i]<6.5:
        ax.vlines(k_llista[i],0,1.1, lw =1, color = 'blue', alpha=0.5)
        if parells[i][0] != 3:
            ax.text(k_llista[i]-0.1, 1.05, r"$K_{"+str(parells[i][0])+","+str(parells[i][1])+"}$", fontsize =9)

ax.plot(dadesT[0][10:]/np.pi ,dadesT[2][10:]/(max(dadesT[2][10:250])), label = 'Placa hexagonal')
ax.plot(dadesT2[0][10:]/np.pi ,dadesT2[2][10:]/(max(dadesT2[2][10:250])), label =' Placa triangular')
ax.legend(loc = 6)
ax.set_xlim(1,6.5)
ax.set_ylim(0,1.1)
ax.set_xlabel(r'$\frac{\hat{k}}{\pi}$',fontsize=22 )
ax.set_ylabel(r'$\frac{E}{E_{max}}$', fontsize=22)

plt.show()


