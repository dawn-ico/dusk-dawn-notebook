import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import glob
from IPython.display import HTML

def update_plot(frame_number, zarray, plot):
    z = np.loadtxt(out_files[frame_number])
    sol = np.loadtxt(sol_files[frame_number])
    for i in range(0,3):
        plot[i].remove()    
    plot[0] = ax1.plot_trisurf(x, y, z, triangles=T-1, cmap=plt.cm.autumn)
    plot[1] = ax2.plot_trisurf(x, y, sol, triangles=T-1, cmap=plt.cm.autumn)
    plot[2] = ax3.plot_trisurf(x, y, z-sol, triangles=T-1, cmap=plt.cm.autumn)
    
zarray=[]

out_files = sorted(glob.glob('out/stepT*.txt'))
sol_files = sorted(glob.glob('out/solT*.txt'))

x,y = np.loadtxt('initP.txt', unpack=True)
T = np.loadtxt('initT.txt')
z = np.loadtxt('initC.txt')
fig = plt.figure(figsize=(12,8), dpi=100)
ax1 = fig.add_subplot(131, projection='3d')  
ax1.set_zlim(0,1)
ax1.set_title('Simulation')
ax2 = fig.add_subplot(132, projection='3d')  
ax2.set_zlim(0,1)
ax2.set_title('Solution')
ax3 = fig.add_subplot(133, projection='3d')  
ax3.set_zlim(-0.1,0.1)
ax3.set_title('Error (Sim - Sol)')
plot = [ax1.plot_trisurf(x, y, np.zeros(np.size(z)), triangles=T-1, cmap=plt.cm.autumn), 
        ax2.plot_trisurf(x, y, np.zeros(np.size(z)), triangles=T-1, cmap=plt.cm.autumn), 
        ax3.plot_trisurf(x, y, np.zeros(np.size(z)), triangles=T-1, cmap=plt.cm.autumn)]