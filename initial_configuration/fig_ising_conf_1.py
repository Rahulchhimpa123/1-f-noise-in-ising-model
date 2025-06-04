import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['text.usetex'] = True

df1 = pd.read_csv('fort.10', delimiter='\\s+', header=None, names=['T11', 'T12'])
df2 = pd.read_csv('fort.11', delimiter='\\s+', header=None, names=['T21', 'T22'])

fig, axd = plt.subplot_mosaic([['F', 'G']], figsize=(3.2,1.6), layout='constrained')


axd['F'].set_ylim(0,128)
axd['F'].set_yticks([0,64,128])
axd['F'].set_xlim(0,128)
axd['F'].set_xticks([0,64,128])
axd['F'].tick_params(axis='x', labelsize=8)
axd['F'].tick_params(axis='y', labelsize=8)


axd['F'].scatter(df1["T11"],df1["T12"],color='C3', marker='+',s=0.25)



axd['G'].set_ylim(0,128)
axd['G'].set_yticks([0,64,128])
axd['G'].set_xlim(0,128)
axd['G'].set_xticks([0,64,128])
axd['G'].tick_params(axis='x', labelsize=8)
axd['G'].tick_params(axis='y', labelsize=8)


axd['G'].scatter(df2["T21"],df2["T22"],color='C5', marker='+',s=0.25)

fig.savefig("fig_ising_conf_1.pdf")
