# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as interpolate



# # For latex font, i guess so 
# plt.rc('text', usetex=True)
plt.rc('font', family='arial')
#Set global matplotlib parameters in script or in /home/$USER/.matplotlib/matplotlibrc
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams.update({'font.size': 20})


xx6,phi6,phi_an6,u6,u_an6=np.loadtxt('Result-1b.plt', delimiter=None, unpack=True)

# xx6,phi6,phi_an6,u6,u_an6=np.loadtxt('7Compact.txt', delimiter=None, unpack=True)

# plt.rc('legend',fontsize=8) # using a size in points
# plt.rc('legend',fontsize='8') # using a named size
# ----------------# 6th order plots------------------------ 





# ----------------# 4th order plots------------------------ 

fg = plt.figure(figsize=(12,6))
sp = plt.subplot(1,2,1)
q1 = plt.plot(xx6, phi6, 'o',markersize=4,markeredgecolor='blue',markerfacecolor='none',label="Numerical solution")
q1 = plt.plot(xx6, phi_an6, '*',markersize=4,markeredgecolor='black',markerfacecolor='none',label="exact solution")
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'$\phi$', fontsize=18)
# plt.legend(loc='best')
# tl = plt.title("Potential-$6^{th}$ order ")

sp = plt.subplot(1,2,2)
q1 = plt.plot(xx6, u6, 'o',markersize=2,markeredgecolor='blue',markerfacecolor='none')
q3 = plt.plot(xx6, u_an6, '*',markersize=2,markeredgecolor='red',markerfacecolor='none')
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'$u$', fontsize=18)

plt.savefig('Nishikawa_2007.pdf',dpi=600)
# plt.show()

