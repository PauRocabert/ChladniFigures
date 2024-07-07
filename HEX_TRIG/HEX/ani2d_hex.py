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
    z[i]=abs(uall[(i*nx):((i+1)*nx)])
zarray=np.transpose(z,(2,1,0))

b = np.linspace(0,1,nx)
a = np.linspace(0,2,ny)

b, a = np.meshgrid(b,a)
x = a
y = b*np.sqrt(3)

num_levels = 100
u_rang = np.linspace(0,np.amax(z), num_levels)

colors = plt.cm.gray(np.linspace(0, 1, num_levels + 1))
colors[:5] = [1, 1, 1, 1] # set to white
colors[6:]=[0,0,0,1]

def update_plot(frame_number, zarray, plot):
    ax.clear()
    ax.contourf(x,y, zarray[:,:,frame_number], num_levels, colors = colors, origin ='lower')

fig, ax = plt.subplots()

plot = [ax.contourf(x,y, zarray[:,:,0], num_levels, colors = colors, origin ='lower')]
ani = animation.FuncAnimation(fig, update_plot, frn, fargs=(zarray, plot), interval=1000/fps)

fn = 'plothex_fig'
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps, dpi =200)

plt.rcParams['animation.html'] = 'html5'
ani
