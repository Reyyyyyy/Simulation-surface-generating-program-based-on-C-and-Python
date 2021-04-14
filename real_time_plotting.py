import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(dpi=128)
ax = Axes3D(fig)
xs = np.arange(-10,10,0.1)
ys = np.arange(-10,10,0.1)
X,Y = np.meshgrid(xs,ys)
zs = []

for n in range(50):
    zs.append((X+0.3*n)**2+(Y-0.3*n)**2+X*Y*n*0.1+n*np.sin(X+Y))

for Z in zs:
    plt.cla()#清除原有图像
    plt.title('Real-time plotting')
    ax.plot_surface(X,Y,Z,cmap=plt.get_cmap('gist_ncar'))
    plt.axis('off')
    plt.pause(0.001)
    

