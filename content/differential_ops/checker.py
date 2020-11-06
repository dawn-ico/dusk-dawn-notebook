import sys
import numpy as np
import matplotlib.pyplot as plt

if (len(sys.argv) != 2):
    print("usage checker gradient|divergence|curl")
    
if (not (str(sys.argv[1]) == "gradient" or str(sys.argv[1]) == "divergence" or str(sys.argv[1]) == "curl")):
    print("usage checker gradient|divergence|curl")
    
T = np.loadtxt('out/meshT.txt')
x,y = np.loadtxt('out/meshP.txt', unpack=True)
    
if (str(sys.argv[1]) == "gradient"):
    f_x = np.loadtxt('out/f_x.txt')
    f_x_Sol = np.loadtxt('out/f_x_Sol.txt')
    f_y = np.loadtxt('out/f_y.txt')
    f_y_Sol = np.loadtxt('out/f_y_Sol.txt')
    fig = plt.figure(figsize=(12,4), dpi=100)
    
    
    ax1 = fig.add_subplot(131)
    tp = ax1.tripcolor(x, y, triangles=T-1, facecolors=f_x)    
    ax1.set_ylim(-2, 2)
    ax1.set_aspect('equal')
    ax1.set_title('gradient_x computed', fontsize=10)    
    plt.colorbar(tp,ax=ax1)
    
    ax2 = fig.add_subplot(132)
    tp = ax2.tripcolor(x, y, triangles=T-1, facecolors=f_x_Sol)
    ax2.set_ylim(-2, 2)
    ax2.set_aspect('equal')
    ax2.set_title('gradient_x analytical Solution', fontsize=10)
    plt.colorbar(tp,ax=ax2)
    
    ax3 = fig.add_subplot(133)
    tp = ax3.tripcolor(x, y, triangles=T-1, facecolors=np.abs(f_x-f_x_Sol))
    ax3.set_ylim(-2, 2)
    ax3.set_aspect('equal')
    plt.colorbar(tp,ax=ax3)
    ax3.set_title('difference (note color scale)', fontsize=10)
    
    plt.subplots_adjust(wspace=0.4)
    plt.show()