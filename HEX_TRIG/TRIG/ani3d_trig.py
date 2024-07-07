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
N = int(len(Lines)/k)
z = np.zeros(shape=(N,m,k))

fps = 10 # frame per sec
frn = N

for n,line in enumerate(Lines):
    for i,u in enumerate(line[:-2].split(' ')):
        z[int(n/k)][i][n%k] = float(u)
zarray=np.transpose(z,(1,2,0))

a = np.linspace(0,1,m)
b = np.linspace(0,0.5,k)

am, bm = np.meshgrid(b,a)
y = am*np.sqrt(3)
x = bm

def update_plot(frame_number):
    global zarray
    ax.clear()
    ax.set_zlim([np.amin(z),np.amax(z)])
    ax.contour3D(x, y, zarray[:,:,frame_number], 50)

fig = plt.figure()
ax = plt.axes(projection='3d')



ax.set_xlim(0,1)
ax.set_ylim(0,1)
ani = animation.FuncAnimation(fig, update_plot, frn, interval=1000/fps)

fn = 'plottrig_3d_'
ani.save(fn+w+'.mp4',writer='ffmpeg',fps=fps,dpi=200)

plt.rcParams['animation.html'] = 'html5'
ani
