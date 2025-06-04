import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['text.usetex'] = True
from matplotlib.ticker import LogLocator


df1 = pd.read_csv('fort.21', delimiter='\\s+', header=None, names=['T11', 'T12'])
df2 = pd.read_csv('fort.22', delimiter='\\s+', header=None, names=['T21', 'T22'])
df3 = pd.read_csv('fort.23', delimiter='\\s+', header=None, names=['T31', 'T32'])
df4 = pd.read_csv('fort.24', delimiter='\\s+', header=None, names=['T41', 'T42'])

fig, axd = plt.subplot_mosaic([['F']], figsize=(3.3,3.3), layout='constrained')#, gridspec_kw={'width_ratios': [1.0],'height_ratios': [1.6,2.7]})

a1 = 2**20-100001



axd['F'].set_xlabel(r'$\tau/N^{z/{\rm D}}$')
axd['F'].set_ylabel(r'$C(\tau,N)/N$')  # Add a y-label to the axes.
#axd['F'].set_yscale('log')
axd['F'].set_xscale('log')
axd['F'].set_ylim(-0.01,0.45)
axd['F'].set_xlim(1e-5,1e4)
axd['F'].set_xticks([1e-5,1e-2,1e1,1e4])
axd['F'].set_yticks([0,0.2,0.4])
axd['F'].text(2e-5, 0.3, '$z = 2.08$', fontsize=9)
#axd['F'].minorticks_off()
axd['F'].xaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=100))
a=-1.0


y=-1.04
s1 = 1.5
N=2**12
axd['F'].plot(s1*(N**(y))*df4["T41"],(N**(a))*(df4["T42"]/a1), label=r'$N=2^{12}$',color='C3', linewidth=1.00, linestyle='-.')

N=2**10
axd['F'].plot(s1*(N**(y))*df3["T31"],(N**(a))*(df3["T32"]/a1), label=r'$N=2^{10}$',color='C2', linewidth=1.00, linestyle='--')#(5,(10,3)))

N=2**8
axd['F'].plot(s1*(N**(y))*df2["T21"],(N**(a))*(df2["T22"]/a1), label=r'$N=2^{8}$',color='C1', linewidth=1.00, linestyle=':')#(0,(5,10)))

N=2**6
axd['F'].plot(s1*(N**(y))*df1["T11"],(N**(a))*(df1["T12"]/a1), label=r'$N=2^{6}$',color='C0', linewidth=1.00, linestyle='-')


d = 1.0
x11 = np.array([2.2e0,6e2])
#axd['F'].plot(x11, (1.1e-1)*x11**(-d), label='',color='black', linewidth=2, linestyle='--')
#axd['F'].text(5e-1, 2e-3, '-1', fontsize=11)

#axd['F'].legend(loc='lower left', frameon=False)  # Add a legend.





left, bottom, width, height = [0.66, 0.36, 0.3, 0.25]
axd['C'] = fig.add_axes([left, bottom, width, height])

axd['C'].xaxis.set_ticks_position('bottom')
axd['C'].xaxis.set_label_position('bottom')
axd['C'].yaxis.set_ticks_position('left')
axd['C'].yaxis.set_label_position('left')


axd['C'].set_xlabel(r'$\tau$',fontsize=7)
axd['C'].set_ylabel(r'$C(\tau, N)/N$',fontsize=7)  # Add a y-label to the axes.
#axd['C'].set_yscale('log')
axd['C'].set_xscale('log')
#axd['C'].set_ylim(1e-4,1e1)
#axd['C'].set_xlim(1e-3,1e5)
axd['C'].set_xlim(1e0,1e6)
axd['C'].set_ylim(-0.01,0.45)
axd['C'].set_xticks([1e0,1e3,1e6])
#axd['C'].set_xticks([0.00001,0.0001,0.001,0.01,0.1,1])
axd['C'].set_yticks([0,0.2,0.4])
#axd['B'].annotate("", xy=(0.5, 5), xytext=(0.05, 0.005),arrowprops=dict(arrowstyle="<-"))
#axd['B'].annotate("", xy=(1e-4, 5), xytext=(1e-4, 1e9),arrowprops=dict(arrowstyle="<-"))
axd['C'].text(1e4, 3.5e-1, r'(b)', fontsize=11)
#axd['C'].xaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=100))
axd['C'].minorticks_on()

a=-2


y=0#-2
s1 = 1

L=2**6
axd['C'].plot(s1*(L**(y))*df4["T41"],(L**(a))*(df4["T42"]/a1), label=r'$N=2^{12}$',color='C3', linewidth=1.00, linestyle='-.')

L=2**5
axd['C'].plot(s1*(L**(y))*df3["T31"],(L**(a))*(df3["T32"]/a1), label=r'$L=2^{5}$',color='C2', linewidth=1.00, linestyle='--')#(5,(10,3)))

L=2**4
axd['C'].plot(s1*(L**(y))*df2["T21"],(L**(a))*(df2["T22"]/a1), label=r'$L=2^{4}$',color='C1', linewidth=1.00, linestyle=':')#(0,(5,10)))

L=2**3
axd['C'].plot(s1*(L**(y))*df1["T11"],(L**(a))*(df1["T12"]/a1), label=r'$L=2^{3}$',color='C0', linewidth=1.00, linestyle='-')


left, bottom, width, height = [0.52, 0.718, 0.4, 0.27]
axd['C'] = fig.add_axes([left, bottom, width, height])

axd['C'].xaxis.set_ticks_position('bottom')
axd['C'].xaxis.set_label_position('bottom')
axd['C'].yaxis.set_ticks_position('left')
axd['C'].yaxis.set_label_position('left')

axd['C'].set_xlabel(r'$\tau$',fontsize=7)
axd['C'].set_ylabel(r'$C(\tau, N)$',fontsize=7)  # Add a y-label to the axes.
#axd['C'].set_yscale('log')
axd['C'].set_xscale('log')
#axd['C'].set_ylim(1e-4,1e1)
#axd['C'].set_xlim(1e-3,1e5)
axd['C'].set_xlim(1e0,1e6)
axd['C'].set_ylim(-100,1.8e3)
axd['C'].set_xticks([1e0,1e3,1e6])
#axd['C'].set_xticks([0.00001,0.0001,0.001,0.01,0.1,1])
#axd['C'].set_yticks([1e-12,1e-7,1e-2])
#axd['B'].annotate("", xy=(0.5, 5), xytext=(0.05, 0.005),arrowprops=dict(arrowstyle="<-"))
#axd['B'].annotate("", xy=(1e-4, 5), xytext=(1e-4, 1e9),arrowprops=dict(arrowstyle="<-"))
axd['C'].text(1e5, 4e1, r'(a)', fontsize=11)
axd['C'].minorticks_on()
#axd['C'].text(2e-2, 1.0e-7, r'$-1.9$', fontsize=7)
a=0#-2


y=0#-2
s1 = 1

L=2**6
axd['C'].plot(s1*(L**(y))*df4["T41"],(L**(a))*(df4["T42"]/a1), label=r'$N=2^{12}$',color='C3', linewidth=1.00, linestyle='-.')

L=2**5
axd['C'].plot(s1*(L**(y))*df3["T31"],(L**(a))*(df3["T32"]/a1), label=r'$N=2^{10}$',color='C2', linewidth=1.00, linestyle='--')#(5,(10,3)))

L=2**4
axd['C'].plot(s1*(L**(y))*df2["T21"],(L**(a))*(df2["T22"]/a1), label=r'$N=2^{8}$',color='C1', linewidth=1.00, linestyle=':')#(0,(5,10)))

L=2**3
axd['C'].plot(s1*(L**(y))*df1["T11"],(L**(a))*(df1["T12"]/a1), label=r'$N=2^{6}$',color='C0', linewidth=1.00, linestyle='-')

axd['C'].legend(loc='upper right', frameon=False, fontsize=8)  # Add a legend.



left, bottom, width, height = [0.24, 0.25, 0.22, 0.22]
axd['C'] = fig.add_axes([left, bottom, width, height])

axd['C'].xaxis.set_ticks_position('bottom')
axd['C'].xaxis.set_label_position('bottom')
axd['C'].yaxis.set_ticks_position('left')
axd['C'].yaxis.set_label_position('left')

df4 = pd.read_csv('fort.59', delimiter='\\s+', header=None, names=['T41', 'T42', 'T43'])

axd['C'].set_xlabel(r'$N$',fontsize=7)
axd['C'].set_ylabel(r'$C(0, N)/N$',fontsize=7)  # Add a y-label to the axes.
#axd['C'].set_yscale('log')
axd['C'].set_xscale('log')
axd['C'].set_xlim(1e1,1e5)
axd['C'].set_ylim(0.2,0.6)
axd['C'].set_xticks([1e1,1e3,1e5])
axd['C'].text(5e3, 3e-1, r'(c)', fontsize=11)
axd['C'].minorticks_on()
a=1
y=0
#axd['C'].plot(s1*(L**(y))*df4["T42"],(df4["T43"]/a1)/((df4["T42"])**(a)),color='C0', linewidth=1.00, linestyle='-')
axd['C'].scatter(s1*(L**(y))*df4["T42"],(df4["T43"]/a1)/((df4["T42"])**(a)),color='C4', marker='^',s=15)

fig.savefig("fig_13.pdf")
#plt.show()

