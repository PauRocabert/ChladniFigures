import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib

w="sqrt5pi"
Zall = np.loadtxt("MNII_PF/QUAD/EX/quad_ex.txt")

fps = 10 # frame per sec
nx=int(len(Zall[0])) #intervalos
frn =  int(len(Zall)/nx) # frame number of the animation

z= np.zeros((frn,nx,nx))

for i in range(frn):
    z[i]=Zall[(i*nx):((i+1)*nx)]

x = np.linspace(-1/2, 1/2, nx)
y = np.linspace(-1/2, 1/2, nx)
zarray=np.transpose(z,(1,2,0))
X, Y = np.meshgrid(x, y)

def update_plot(frame_number):
    global zarray
    ax.clear()
    ax.set_zlim([np.amin(Zall),np.amax(Zall)])
    ax.contour3D(x, y, zarray[:,:,frame_number], 50)

fig = plt.figure()
ax = plt.axes(projection='3d')

plot = [ax.contour3D(x, y, zarray[:,:,0], 50)]

ani = animation.FuncAnimation(fig, update_plot, frn, interval=1000/fps)

fn = 'plotquad_3d_'+w
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps,dpi=200)

plt.rcParams['animation.html'] = 'html5'
ani
