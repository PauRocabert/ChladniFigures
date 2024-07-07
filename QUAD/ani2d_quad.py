import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib

Zall = np.loadtxt("MNII_PF/QUAD/EX/quad_ex.txt")

w="sqrt5pi_ex"

fps = 10 # frame per sec
nx=int(len(Zall[0])) #intervalos
frn = int(len(Zall)/nx) # frame number of the animation

z= np.zeros((frn,nx,nx))

for i in range(frn):
    q=Zall[(i*nx):((i+1)*nx)]
    z[i]=abs(q)

x = np.linspace(-0.5, 0.5, nx)
y = np.linspace(-0.5, 0.5, nx)
zarray=np.transpose(z,(1,2,0))
X, Y = np.meshgrid(x, y)

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

fn = 'plotquad_fig_'+w
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps, dpi =200)

plt.rcParams['animation.html'] = 'html5'
ani
