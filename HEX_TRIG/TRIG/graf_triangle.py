import numpy as np 
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt 

freq = '(11)'
file = open(f'MNII_PF/HEX_TRIG/TRIG/trig.txt')

Lines = file.readlines()

m =len(Lines[0][:-2].split(' '))
k = int((m+1)/2)
N = int(len(Lines)/k)
u = np.zeros(shape=(N,m,k))

for n,line in enumerate(Lines):
    for i,z in enumerate(line[:-2].split(' ')):
        u[int(n/k)][i][n%k] = float(z)
        
uT = np.transpose(u)

a = np.linspace(0,1,m)
b = np.linspace(0,0.5,k)


am, bm = np.meshgrid(b,a)
y = am*np.sqrt(3)
x = bm

ax = plt.axes(projection='3d')
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.contour3D(x,y,u[-1], 50, cmap='viridis')
plt.show()

#chladni fig##################################################################################################

ur=abs(u[-1])
umax = np.max(ur)
umin = np.min(ur)
num_levels = 100
u_rang = np.linspace(0,umax, num_levels)

colors = plt.cm.gray(np.linspace(0, 1, num_levels + 1))
colors[:5] = [1, 1, 1, 1] # set to white
colors[6:] = [0, 0, 0, 1]


fig, ax = plt.subplots()
ax.contourf(x,y, ur, num_levels, colors = colors, origin ='lower')
ax.set_xlabel(r'$\hat{x}$', fontsize=16)
ax.set_ylabel(r'$\hat{y}$', fontsize=16)
plt.show()