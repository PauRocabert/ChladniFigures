import numpy as np 
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt 

uall = np.loadtxt('MNII_PF/HEX_TRIG/HEX/hex.txt')

ny =len(uall[0])
nx = int(ny/2)
t = int(len(uall)/nx)
u = np.zeros(shape=(t,nx,ny))

for i in range(t):
    u[i]=uall[(i*nx):((i+1)*nx)]
u=np.transpose(u,(0,2,1))

b = np.linspace(0,1,nx)
a = np.linspace(0,2,ny)

b, a = np.meshgrid(b,a)
x = a
y = b*np.sqrt(3)


ax = plt.axes(projection='3d')
ax.contour3D(x,y,u[-1], 50, cmap='viridis')
ax.set_xlabel(r'$\hat{x}$')
ax.set_ylabel(r'$\hat{y}$')

plt.show()

#chladni fig##################################################################################################
ur=abs(u[-1])
umax = np.max(ur)
umin = np.min(ur)
num_levels = 100
u_rang = np.linspace(0,umax, num_levels)

colors = plt.cm.gray(np.linspace(0, 1, num_levels + 1))
colors[:3] = [1, 1, 1, 1] # set to white
colors[6:] = [0, 0, 0, 1]


fig, ax = plt.subplots(figsize=(8,8))
ax.contourf(x,y, ur, num_levels, colors = colors, origin ='lower')
ax.set_xlabel(r'$\hat{x}$', fontsize=22)
ax.set_ylabel(r'$\hat{y}$', fontsize=22)
plt.show()