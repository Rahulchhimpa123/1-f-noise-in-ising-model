import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['text.usetex'] = True
from matplotlib.ticker import LogLocator


df2 = pd.read_csv('fort.21', delimiter='\\s+', header=None, names=['T21', 'T22'])
df3 = pd.read_csv('fort.22', delimiter='\\s+', header=None, names=['T31', 'T32'])
df4 = pd.read_csv('fort.23', delimiter='\\s+', header=None, names=['T41', 'T42'])
df5 = pd.read_csv('fort.24', delimiter='\\s+', header=None, names=['T51', 'T52'])
df6 = pd.read_csv('fort.25', delimiter='\\s+', header=None, names=['T61', 'T62'])

fig, axd = plt.subplot_mosaic([['F']], figsize=(3.3,3.3), layout='constrained')



axd['F'].set_xlabel(r'$fN^{z/{\rm D}}$')
#axd['F'].set_xlabel(r'$f$')
axd['F'].set_ylabel(r'$S(f,N)/N^{1+{z/{\rm D}}}$')  # Add a y-label to the axes.
axd['F'].set_yscale('log')
axd['F'].set_xscale('log')
axd['F'].set_ylim(1e-5,1.4e-1)
axd['F'].set_xlim(1e-5,1e4)
axd['F'].set_xticks([1e-5,1e-2,1e1,1e4])
#axd['F'].set_yticks([1e-2,1e3,1e8,1e13])
axd['F'].text(2e-5, 1e-2, '$z = 2.08$', fontsize=9)
axd['F'].minorticks_on()
axd['F'].xaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=100))

a=-2.04#-3#-5.8#-6#3.35 # tau b = 2
a1 = 2.67

y=1.04 #lambda or D

s1 = 1.5e0

N=2**12
axd['F'].plot(s1*(N**(y))*df6["T61"],(N**(a))*df6["T62"], label=r'$N=2^{12}$',color='C4', linewidth=1.00, linestyle='-')
N=2**10
axd['F'].plot(s1*(N**(y))*df5["T51"],(N**(a))*df5["T52"], label=r'$N=2^{10}$',color='C3', linewidth=1.00, linestyle=(0,(5,10)))
N=2**8
axd['F'].plot(s1*(N**(y))*df4["T41"],(N**(a))*df4["T42"], label=r'$N=2^{8}$',color='C2', linewidth=1.00, linestyle='-.')

N=2**6
axd['F'].plot(s1*(N**(y))*df3["T31"],(N**(a))*df3["T32"], label=r'$N=2^{6}$',color='C1', linewidth=1.00, linestyle='--')

N=2**4
axd['F'].plot(s1*(N**(y))*df2["T21"],(N**(a))*df2["T22"], label=r'$N=2^{4}$',color='C0', linewidth=1.00, linestyle=':')



d = 1.0
x11 = np.array([1e0,1e3])
axd['F'].plot(x11, (4.5e-2)*x11**(-d), label='',color='black', linewidth=2, linestyle='--')
axd['F'].text(7e-1, 4e-3, '-1', fontsize=11)

#axd['F'].legend()  # Add a legend.
axd['F'].legend(loc='upper right',frameon=False,  fontsize=9)

left, bottom, width, height = [0.334, 0.25, 0.38, 0.37]
axd['C'] = fig.add_axes([left, bottom, width, height])

axd['C'].xaxis.set_ticks_position('bottom')
axd['C'].xaxis.set_label_position('bottom')
axd['C'].yaxis.set_ticks_position('left')
axd['C'].yaxis.set_label_position('left')

axd['C'].set_xlabel(r'$f$',fontsize=7)
axd['C'].set_ylabel(r'$S(f, N)$',fontsize=7)  # Add a y-label to the axes.
axd['C'].set_yscale('log')
axd['C'].set_xscale('log')
axd['C'].set_ylim(1e0,3e6)
axd['C'].set_xlim(1e-6,1e0)
#axd['C'].set_xlim(1e0,1e4)
#axd['C'].set_yticks([1e0,1e3,1e6])
axd['C'].set_xticks([1e-6,1e-3,1e0])

axd['C'].annotate("", xy=(5e-6, 5e0), xytext=(5e-6, 3e6),
            arrowprops=dict(arrowstyle="<-"))


a=0

y=0#2.0 #lambda or D

s1 = 1e0

L=2**6
axd['C'].plot(s1*(L**(y))*df6["T61"],(L**(a))*df6["T62"], label=r'$N=2^{12}$',color='C4', linewidth=1.00, linestyle='-')

L=2**5
axd['C'].plot(s1*(L**(y))*df5["T51"],(L**(a))*df5["T52"], label=r'$N=2^{10}$',color='C3', linewidth=1.00, linestyle=(0,(5,10)))
L=2**4
axd['C'].plot(s1*(L**(y))*df4["T41"],(L**(a))*df4["T42"], label=r'$N=2^{8}$',color='C2', linewidth=1.00, linestyle='-.')

L=2**3
axd['C'].plot(s1*(L**(y))*df3["T31"],(L**(a))*df3["T32"], label=r'$N=2^{6}$',color='C1', linewidth=1.00, linestyle='--')

L=2**2
axd['C'].plot(s1*(L**(y))*df2["T21"],(L**(a))*df2["T22"], label=r'$N=2^{4}$',color='C0', linewidth=1.00, linestyle=':')

fig.savefig("fig_ising_2d_accepted_count_11.pdf")
#plt.show()

