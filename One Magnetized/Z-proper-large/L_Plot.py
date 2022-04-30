from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


# For latex font, i guess so 
plt.rc('text', usetex=True)
plt.rc('font', family='arial')
plt.rc('legend',fontsize=18) # using a size in points


#Set global matplotlib parameters in script or in /home/$USER/.matplotlib/matplotlibrc
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.size'] = 8
# plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
# plt.rcParams['ytick.minor.size'] = 3
plt.rcParams.update({'font.size': 22})


p=np.loadtxt('Potential.txt', delimiter=None, unpack=True)
u=np.loadtxt('U_x.txt', delimiter=None, unpack=True)
v=np.loadtxt('U_y.txt', delimiter=None, unpack=True)
ntweno,res_pweno,res_uweno,res_vweno=np.loadtxt('residual.txt', delimiter=None, unpack=True)

plt.plot(ntweno, np.log10(res_uweno/res_uweno[0]), '-^',label=r"$\tilde{u}$-WENO-5Z",color='k',markersize=2,markeredgecolor='green',markerfacecolor='none')
plt.plot(ntweno, np.log10(res_pweno/res_pweno[0]), '-o',label=r"$\tilde{\phi}$-WENO-5Z",color='k',markersize=2,markeredgecolor='red',markerfacecolor='none')
plt.plot(ntweno, np.log10(res_vweno/res_vweno[0]), '-o',label=r"$\tilde{v}$-WENO-5Z",color='k',markersize=2,markeredgecolor='blue',markerfacecolor='none')
plt.ylabel(r'Log10(ResA)',size=22)
plt.xlabel(r'Number of iterations',size=22)
plt.legend(loc='best',fontsize=16,frameon=False)
fig2 = plt.gcf()
fig2.set_size_inches(w=6,h=6)
plt.savefig('residual_potential.pdf',dpi=600,bbox_inches='tight', pad_inches = 0)
plt.show()



print(np.min(p))
print(np.max(p))
print(np.min(u))
print(np.max(u))

print(np.min(v))
print(np.max(v))




m=2
g=48*m	
k=48*m
x = np.linspace(0.0,200.0,g)
y = np.linspace(0.0,100.0,k)

fg = plt.figure(figsize=(8,4))
plt.contourf(x,y,p,cmap='jet')
# CS=plt.contour(x,y,p,10,colors='k', linestyles='dashed')

# CS=plt.contour(x,y,p,10,cmap='jet', linestyles='dashed')
# plt.clabel(CS,fontsize=18, inline=True,fmt = '%1.2f')

# plt.colorbar()
plt.xlabel(r'$\tilde{x}$', fontsize=22)
plt.ylabel(r'$\tilde{y}$', fontsize=22)
fg.suptitle(r'$\tilde{\phi}_{max}=$ %s \quad $\tilde{\phi}_{min}=$ %s' %((np.max(p)), (np.min(p))), fontsize=22, fontweight='bold')
plt.savefig('Potential_large.pdf',dpi=600,bbox_inches='tight', pad_inches = 0)
plt.show()


# # X, Y = np.meshgrid(x, y)
# # fig = plt.figure(figsize=(8,8))
# # ax = fig.gca(projection='3d')
# # cset = ax.plot_surface(X, Y, v, cmap=cm.jet,linewidth=0, antialiased=False)
# # # ax.clabel(cset, fontsize=9, inline=1)
# # # plt.axis('off')
# # plt.xlabel(r'$x$', fontsize=18)
# # plt.ylabel(r'$y$', fontsize=18)
# # # fig.colorbar(cset, shrink=0.5, aspect=5)
# # for angle in range(0, 360):
# #     ax.view_init(elev=30., azim=66.0)
# #     plt.draw()
# # plt.savefig('3D_disc_good.pdf',dpi=600,bbox_inches='tight', pad_inches = 0.05)
# # plt.show()


# fg = plt.figure(figsize=(10,6.5))
fg = plt.figure(figsize=(8,4))
CS=plt.contour(x,y,v,10,cmap='jet')
plt.colorbar()
plt.xlabel(r'$\tilde{x}$', fontsize=22)
plt.ylabel(r'$\tilde{y}$', fontsize=22)
plt.savefig('V_NEW_WENO_1000.pdf',dpi=600,bbox_inches='tight', pad_inches = 0)
plt.show()

fg = plt.figure(figsize=(8,4))
CS=plt.contour(x,y,u,10,cmap='jet')
plt.colorbar()
plt.xlabel(r'$\tilde{x}$', fontsize=22)
plt.ylabel(r'$\tilde{y}$', fontsize=22)
plt.savefig('U_NEW_WENO_1000.pdf',dpi=600,bbox_inches='tight', pad_inches = 0)
plt.show()


fg = plt.figure(figsize=(8,4))
plt.streamplot(x, y, u, v, density = 1.0, color = 'k', arrowsize = 1, arrowstyle = '->', minlength = 0.55, linewidth = 0.65, transform = None, cmap=plt.cm.autumn)
plt.xlabel(r'$\tilde{x}$', fontsize=22)
plt.ylabel(r'$\tilde{y}$', fontsize=22)
plt.xlim(0.0,200.0)
plt.ylim(0.0,100.0)
plt.savefig('Electron_largest.pdf',dpi=600,bbox_inches='tight', pad_inches = 0)
plt.show()







