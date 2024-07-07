import numpy as np 
import matplotlib.pyplot as plt 

u = np.loadtxt('MNII_PF/RAD/rad.txt')
theta = np.linspace(0,2*np.pi, 1000)
r = np.linspace(0,1, len(u[0]))

r,thetaM = np.meshgrid(r, theta)
x = np.cos(thetaM)*r
y = np.sin(thetaM)*r

print(len(x),len(y))

uM, thetaM2 = np.meshgrid(u[-1], theta)

ax = plt.axes(projection='3d')
ax.contour3D(x,y,uM, 50, cmap='viridis')
plt.show()

#chladni fig##################################################################################################
ur=abs(uM)
umax = np.max(ur)

num_levels = 100
u_rang = np.linspace(0,umax, num_levels)

colors = plt.cm.gray(np.linspace(0, 1, num_levels + 1))
colors[:5] = [1, 1, 1, 1] # set to white
colors[6:]=[0,0,0,1]

fig, ax = plt.subplots(figsize=(7,7))
ax.contourf(x,y, ur, num_levels, colors = colors, origin ='lower')
ax.set_xlabel(r'$\hat{x}$', fontsize=22)
ax.set_ylabel(r'$\hat{y}$', fontsize=22)
plt.show()