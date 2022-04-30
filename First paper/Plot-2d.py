# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as interpolate


# # For latex font, i guess so 
plt.rc('text', usetex=True)
plt.rc('font', family='arial')
#Set global matplotlib parameters in script or in /home/$USER/.matplotlib/matplotlibrc
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams.update({'font.size': 20})


# Chose your solution file from the folder. Don't be dumb and complain later.

x,y,p,u,v,p_an=np.loadtxt('Result-6a.plt', delimiter=None, unpack=True,skiprows=3)

# We are plotting j and then i. So, if your mesh is NX by NY and then plot NY by NX
# Fortran is column major and it faster that way. ok?
g=16*2
k=16*2



x = x.reshape(g,k)
y = y.reshape(g,k)
p = p.reshape(g,k)
p_an = p_an.reshape(g,k)
u = u.reshape(g,k)
v = v.reshape(g,k)

fg = plt.figure(figsize=(12,6))
sp = plt.subplot(1,2,1)
plt.contourf(x,y,p,cmap='jet')
plt.colorbar()
plt.xlabel(r'x', fontsize=16)
plt.ylabel(r'y', fontsize=16)
# fig = plt.gcf()
# plt.savefig('3_hall.pdf',dpi=600)

sp = plt.subplot(1,2,2)
# fig = plt.figure(figsize=(12,6), dpi=500)
# plt.xlim(0.0,1.0)
# plt.ylim(0.0,1.0)
# plt.streamplot(x, y, u, v, density=3, arrowsize=1, arrowstyle='->',color='b')
plt.contourf(x,y,p_an,cmap='jet')
plt.colorbar()
plt.savefig('5a.pdf',dpi=600)
# plt.show()

# X, Y = np.meshgrid(x, y)
# fig = plt.figure(figsize=(8,8))
# ax = fig.gca(projection='3d')
# cset = ax.plot_surface(x, y, p, cmap=cm.jet,linewidth=0, antialiased=False)
# # ax.clabel(cset, fontsize=9, inline=1)
# # plt.axis('off')
# plt.xlabel(r'$x$', fontsize=18)
# plt.ylabel(r'$y$', fontsize=18)
# fig.colorbar(cset, shrink=0.5, aspect=5)
# # plt.savefig('3D_exact.pdf',dpi=600)
# plt.show()