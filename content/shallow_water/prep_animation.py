import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import glob
from IPython.display import HTML

def update_plot(frame_number, zarray, plot):
    z = np.loadtxt(out_files[frame_number])    
    plot[0].remove()    
    plot[0] = ax1.plot_trisurf(x, y, z, triangles=T-1, cmap=plt.cm.winter)    
    
zarray=[]

out_files = sorted(glob.glob('out/stepH*.txt'))

x,y = np.loadtxt('out/initP.txt', unpack=True)
T = np.loadtxt('out/initT.txt')
z = np.loadtxt('out/initC.txt')
fig = plt.figure(figsize=(12,8), dpi=100)
ax1 = fig.add_subplot(111, projection='3d')  
ax1.set_zlim(1.5,2.5)
ax1.set_title('Shallow Water')
plot = [ax1.plot_trisurf(x, y, np.zeros(np.size(z)), triangles=T-1, cmap=plt.cm.autumn)]