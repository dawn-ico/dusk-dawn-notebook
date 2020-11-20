import sys
import os
import numpy as np
import matplotlib.pyplot as plt

if (len(sys.argv) != 2):
    print("usage checker gradient|divergence|curl")
    sys.exit()    
    
if (not (str(sys.argv[1]) == "gradient" or str(sys.argv[1]) == "divergence" or str(sys.argv[1]) == "curl" or str(sys.argv[1]) == "laplacian")):
    print("usage checker gradient|divergence|curl|laplacian")
    sys.exit()
    
def file_exists(PATH):
    return os.path.isfile(PATH) and os.access(PATH, os.R_OK)
    
    
if not file_exists('out/meshT.txt') or not file_exists('out/meshP.txt'):
    print("please make sure that you ran the runner before executing the checker")
    sys.exit()
    
T = np.loadtxt('out/meshT.txt')
x,y = np.loadtxt('out/meshP.txt', unpack=True)

if (str(sys.argv[1]) == "laplacian"):
    if not file_exists('out/uv_nabla2.txt'):
        print("please make sure that you ran ./runner laplacian before checking the laplacian")
        sys.exit()
    
    laplacian_uv = np.loadtxt('out/uv_nabla2.txt')    
    #curl_uv_Sol = np.loadtxt('out/uv_curl_Sol.txt')    
    fig = plt.figure(figsize=(12,4), dpi=100)    
    
    ax1 = fig.add_subplot(131)
    tp = ax1.tripcolor(x, y, triangles=T-1, facecolors=laplacian_uv)    
    ax1.set_ylim(-2, 2)
    ax1.set_aspect('equal')
    ax1.set_title('laplacian', fontsize=10)    
    plt.colorbar(tp,ax=ax1)
    
    #ax2 = fig.add_subplot(132)
    #tp = ax2.tripcolor(x, y, triangles=T-1, facecolors=curl_uv_Sol)
    #ax2.set_ylim(-2, 2)
    #ax2.set_aspect('equal')
    #ax2.set_title('divergence solution', fontsize=10)
    #plt.colorbar(tp,ax=ax2)
    
    #ax3 = fig.add_subplot(133)
    #tp = ax3.tripcolor(x, y, triangles=T-1, facecolors=np.abs(curl_uv-curl_uv_Sol))
    #ax3.set_ylim(-2, 2)
    #ax3.set_aspect('equal')
    #plt.colorbar(tp,ax=ax3)
    #ax3.set_title('difference (note color scale)', fontsize=10)

if (str(sys.argv[1]) == "curl"):
    if not file_exists('out/uv_curl.txt') or not file_exists('out/uv_curl_Sol.txt'):
        print("please make sure that you ran ./runner curl before checking the curl")
        sys.exit()
    
    curl_uv = np.loadtxt('out/uv_curl.txt')
    curl_uv_Sol = np.loadtxt('out/uv_curl_Sol.txt')    
    fig = plt.figure(figsize=(12,4), dpi=100)    
    
    ax1 = fig.add_subplot(131)
    tp = ax1.tripcolor(x, y, triangles=T-1, facecolors=curl_uv)    
    ax1.set_ylim(-2, 2)
    ax1.set_aspect('equal')
    ax1.set_title('divergence', fontsize=10)    
    plt.colorbar(tp,ax=ax1)
    
    ax2 = fig.add_subplot(132)
    tp = ax2.tripcolor(x, y, triangles=T-1, facecolors=curl_uv_Sol)
    ax2.set_ylim(-2, 2)
    ax2.set_aspect('equal')
    ax2.set_title('divergence solution', fontsize=10)
    plt.colorbar(tp,ax=ax2)
    
    ax3 = fig.add_subplot(133)
    tp = ax3.tripcolor(x, y, triangles=T-1, facecolors=np.abs(curl_uv-curl_uv_Sol))
    ax3.set_ylim(-2, 2)
    ax3.set_aspect('equal')
    plt.colorbar(tp,ax=ax3)
    ax3.set_title('difference (note color scale)', fontsize=10)

if (str(sys.argv[1]) == "divergence"):
    if not file_exists('out/uv_div.txt') or not file_exists('out/uv_div_Sol.txt'):
        print("please make sure that you ran ./runner divergence before checking the divergence")
        sys.exit()
    
    div_uv = np.loadtxt('out/uv_div.txt')
    div_uv_Sol = np.loadtxt('out/uv_div_Sol.txt')    
    fig = plt.figure(figsize=(12,4), dpi=100)    
    
    ax1 = fig.add_subplot(131)
    tp = ax1.tripcolor(x, y, triangles=T-1, facecolors=div_uv)    
    ax1.set_ylim(-2, 2)
    ax1.set_aspect('equal')
    ax1.set_title('divergence', fontsize=10)    
    plt.colorbar(tp,ax=ax1)
    
    ax2 = fig.add_subplot(132)
    tp = ax2.tripcolor(x, y, triangles=T-1, facecolors=div_uv_Sol)
    ax2.set_ylim(-2, 2)
    ax2.set_aspect('equal')
    ax2.set_title('divergence solution', fontsize=10)
    plt.colorbar(tp,ax=ax2)
    
    ax3 = fig.add_subplot(133)
    tp = ax3.tripcolor(x, y, triangles=T-1, facecolors=np.abs(div_uv-div_uv_Sol))
    ax3.set_ylim(-2, 2)
    ax3.set_aspect('equal')
    plt.colorbar(tp,ax=ax3)
    ax3.set_title('difference (note color scale)', fontsize=10)
    
if (str(sys.argv[1]) == "gradient"):
    if not file_exists('out/f_x.txt') or not file_exists('out/f_y.txt') or not file_exists('out/f_x_Sol.txt') or not file_exists('out/f_y_Sol.txt'):
        print("please make sure that you ran ./runner gradient before checking the gradient")
        sys.exit()
        
    f_x = np.loadtxt('out/f_x.txt')
    f_x_Sol = np.loadtxt('out/f_x_Sol.txt')
    f_y = np.loadtxt('out/f_y.txt')
    f_y_Sol = np.loadtxt('out/f_y_Sol.txt')
    fig = plt.figure(figsize=(12,8), dpi=100)    
    
    ax1 = fig.add_subplot(231)
    tp = ax1.tripcolor(x, y, triangles=T-1, facecolors=f_x)    
    ax1.set_ylim(-2, 2)
    ax1.set_aspect('equal')
    ax1.set_title('gradient_x computed', fontsize=10)    
    plt.colorbar(tp,ax=ax1)
    
    ax2 = fig.add_subplot(232)
    tp = ax2.tripcolor(x, y, triangles=T-1, facecolors=f_x_Sol)
    ax2.set_ylim(-2, 2)
    ax2.set_aspect('equal')
    ax2.set_title('gradient_x analytical solution', fontsize=10)
    plt.colorbar(tp,ax=ax2)
    
    ax3 = fig.add_subplot(233)
    tp = ax3.tripcolor(x, y, triangles=T-1, facecolors=np.abs(f_x-f_x_Sol))
    ax3.set_ylim(-2, 2)
    ax3.set_aspect('equal')
    plt.colorbar(tp,ax=ax3)
    ax3.set_title('difference (note color scale)', fontsize=10)
    
    ax4 = fig.add_subplot(234)
    tp = ax4.tripcolor(x, y, triangles=T-1, facecolors=f_y)    
    ax4.set_ylim(-2, 2)
    ax4.set_aspect('equal')
    ax4.set_title('gradient_y computed', fontsize=10)    
    plt.colorbar(tp,ax=ax4)
    
    ax5 = fig.add_subplot(235)
    tp = ax5.tripcolor(x, y, triangles=T-1, facecolors=f_y_Sol)
    ax5.set_ylim(-2, 2)
    ax5.set_aspect('equal')
    ax5.set_title('gradient_y analytical solution', fontsize=10)
    plt.colorbar(tp,ax=ax5)
    
    ax6 = fig.add_subplot(236)
    tp = ax6.tripcolor(x, y, triangles=T-1, facecolors=np.abs(f_y-f_y_Sol))
    ax6.set_ylim(-2, 2)
    ax6.set_aspect('equal')
    plt.colorbar(tp,ax=ax6)
    ax6.set_title('difference (note color scale)', fontsize=10)
    
    plt.subplots_adjust(wspace=0.4)
    plt.show()