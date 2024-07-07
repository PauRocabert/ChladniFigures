import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib

u = np.loadtxt('MNII_PF/RAD/rad.txt', delimiter=' ')

fps = 10 # frame per sec
frn = int(len(u)) # frame number of the animation
nr=len(u[0])
nt = 100

w=""
z= np.zeros((frn,nt,nr))

theta = np.linspace(0,2*np.pi, nt)
r = np.linspace(0,1, nr)
r,thetaM = np.meshgrid(r, theta)
x = np.cos(thetaM)*r
y = np.sin(thetaM)*r

z=np.zeros((frn,nt,nr))

for i in range(0,frn):
    uM, thetaM2 = np.meshgrid(u[i], theta)
    z[i]=uM

zarray=np.transpose(z,(1,2,0))

def update_plot(frame_number):
    global zarray
    ax.clear()
    ax.set_zlim([np.amin(z),np.amax(z)])
    ax.contour3D(x, y, zarray[:,:,frame_number], 50)

fig = plt.figure()
ax = plt.axes(projection='3d')

plot = [ax.contour3D(x, y, zarray[:,:,0], 50)]
ani = animation.FuncAnimation(fig, update_plot, frn, interval=1000/fps)

fn = 'plotrad_3d_'+w
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps, dpi=200)

plt.rcParams['animation.html'] = 'html5'
ani
