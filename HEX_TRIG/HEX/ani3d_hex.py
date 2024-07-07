import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib

uall = np.loadtxt('MNII_PF/HEX_TRIG/HEX/hex.txt')

ny =len(uall[0])
nx = int(ny/2)

fps = 10 # frame per sec
frn = int(len(uall)/nx)

z = np.zeros(shape=(frn,nx,ny))

for i in range(frn):
    z[i]=uall[(i*nx):((i+1)*nx)]
zarray=np.transpose(z,(2,1,0))

b = np.linspace(0,1,nx)
a = np.linspace(0,2,ny)

b, a = np.meshgrid(b,a)
x = a
y = b*np.sqrt(3)

def update_plot(frame_number):
    global zarray
    ax.clear()
    ax.axis('off')
    ax.set_zlim([np.amin(uall),np.amax(uall)])
    ax.contour3D(x, y, zarray[:,:,frame_number], 50)

fig = plt.figure()
ax = plt.axes(projection='3d')


ax.set_zlim([np.amin(uall),np.amax(uall)])

ani = animation.FuncAnimation(fig, update_plot, frames = frn, interval=1000/fps)


fn = 'plothex_3d'
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps, dpi=200)
plt.rcParams['animation.html'] = 'html5'
ani
