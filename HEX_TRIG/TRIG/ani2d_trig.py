import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib

file = open(f'MNII_PF/HEX_TRIG/TRIG/trig.txt')
Lines = file.readlines()
w="11"

m =len(Lines[0][:-2].split(' '))
k = int((m+1)/2)
fps = 10 # frame per sec
frn = int(len(Lines)/k)
u = np.zeros(shape=(frn,m,k))

for n,line in enumerate(Lines):
    for i,z in enumerate(line[:-2].split(' ')):
        u[int(n/k)][i][n%k] = abs(float(z))

zarray=np.transpose(u,(1,2,0))

a = np.linspace(0,1,m)
b = np.linspace(0,0.5,k)


am, bm = np.meshgrid(b,a)
y = am*np.sqrt(3)
x = bm

ax = plt.axes(projection='3d')
ax.set_xlim(0,1)
ax.set_ylim(0,1)

num_levels = 100
u_rang = np.linspace(0,np.amax(zarray), num_levels)

colors = plt.cm.gray(np.linspace(0, 1, num_levels + 1))
colors[:5] = [1, 1, 1, 1] # set to white
colors[6:]=[0,0,0,1]

def update_plot(frame_number, zarray, plot):
    ax.clear()
    ax.contourf(x,y, zarray[:,:,frame_number], num_levels, colors = colors, origin ='lower')

fig, ax = plt.subplots()
ax.axis('off')
plot = [ax.contourf(x,y, zarray[:,:,0], num_levels, colors = colors, origin ='lower')]
ani = animation.FuncAnimation(fig, update_plot, frn, fargs=(zarray, plot), interval=1000/fps)

fn = 'plottrig_fig'
ani.save(fn+"_"+w+'.mp4',writer='ffmpeg',fps=fps,dpi=200)

plt.rcParams['animation.html'] = 'html5'
ani
