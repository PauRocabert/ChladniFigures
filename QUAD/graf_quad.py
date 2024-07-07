import numpy as np 
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt 

uall = np.loadtxt("MNII_PF/QUAD/EX/quad_ex.txt")
nx=int(len(uall[0]))

t=int(len(uall)/nx)
u= np.zeros((t,nx,nx))

for i in range(t):
    u[i]=uall[(i*nx):((i+1)*nx)]

#PLOT 3D#######################################################################################################
x = np.linspace(-0.5, 0.5, nx)
y = np.linspace(-0.5, 0.5, nx)
X, Y = np.meshgrid(x, y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.contour3D(X,Y,u[-1], 50, cmap='viridis')
plt.show()

#chladni points##################################################################################################

uq=abs(u[-1])
umax = np.max(u[-1])

num_levels = 100
u_rang = np.linspace(0,umax, num_levels)


colors = plt.cm.gray(np.linspace(0, 1, num_levels + 1))
colors[:5] = [1, 1, 1, 1] # set to white
colors[6:]=[0,0,0,1]


fig, ax = plt.subplots(figsize=(7,7))

ax.contourf(x,y, uq, num_levels, colors = colors, origin ='lower')
ax.set_xlabel(r'$\hat{x}$', fontsize=22)
ax.set_ylabel(r'$\hat{y}$', fontsize=22)
fig.savefig('quadrat.png')
plt.show()
                
