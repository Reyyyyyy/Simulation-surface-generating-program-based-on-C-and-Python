import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#hole
fig_1 = plt.figure(figsize=(8,5),dpi=128)
ax1 = Axes3D(fig_1)
xs_1 = np.arange(-5,5,0.01)
ys_1 = np.arange(-5,5,0.01)
xs1_grid,ys1_grid = np.meshgrid(xs_1,ys_1)
zs_1 = xs1_grid**2 + ys1_grid**2
cmap_1 = plt.get_cmap('rainbow')
ax1.plot_surface(xs1_grid,ys1_grid,zs_1,cmap=cmap_1)
plt.axis('off')

#cap
fig_2 = plt.figure(figsize=(8,5),dpi=128)
ax2 = Axes3D(fig_2)
xs_2 = np.arange(-5,5,0.01)
ys_2 = np.arange(-5,5,0.01)
xs2_grid,ys2_grid = np.meshgrid(xs_2,ys_2)
zs_2 = np.sin(np.sqrt(xs2_grid**2+ys2_grid**2))
cmap_2 = plt.get_cmap('gist_ncar')
ax2.plot_surface(xs2_grid,ys2_grid,zs_2,cmap=cmap_2)
plt.axis('off')

plt.show()
