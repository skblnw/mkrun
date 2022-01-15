#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

params = {
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 16,
    'axes.labelsize': 24,
    'axes.linewidth': 2,
    'axes.xmargin': 0,
    'lines.linewidth' : 1,
    'legend.fontsize': 16,
    'xtick.labelsize': 16,
    'xtick.major.size': 2,
    'xtick.major.width': 2,
    'ytick.labelsize': 16,
    'ytick.major.size': 2,
    'ytick.major.width': 2,
    'text.usetex': False,
    'figure.figsize': [12, 8],
    'figure.dpi': 100
   }
matplotlib.rcParams.update(params)
cmap = matplotlib.cm.get_cmap('Pastel2')

dayback = 10

Jobs={}
UserCount={}
UserRTime={}
UserQTime={}
Colors={}
with open('last'+str(dayback)) as f:
    next(f)
    for line in f:
        (jobid, qname, qtime, rtime, user) = line.split()
        Jobs[int(jobid)] = {'qname': qname, 'qtime': float(qtime), 'rtime': float(rtime), 'user': user}
        UserCount[user] = {'total': 0, 'md': 0, 'fep': 0}
        UserRTime[user] = {'total': 0, 'md': 0, 'fep': 0}
        UserQTime[user] = {'total': 0, 'md': 0, 'fep': 0}
        if user == "flei":
            Colors[user] = 'grey'
        else:
            Colors[user] = 'grey'


for jobid, item in Jobs.items():
    if Jobs[jobid]['qname'] == "md":
        UserCount[Jobs[jobid]['user']][Jobs[jobid]['qname']] += 1
        UserRTime[Jobs[jobid]['user']][Jobs[jobid]['qname']] += Jobs[jobid]['rtime']
        UserQTime[Jobs[jobid]['user']][Jobs[jobid]['qname']] += Jobs[jobid]['qtime']
    else:
        UserCount[Jobs[jobid]['user']][Jobs[jobid]['qname']] += 1
        UserRTime[Jobs[jobid]['user']][Jobs[jobid]['qname']] += Jobs[jobid]['rtime']
        UserQTime[Jobs[jobid]['user']][Jobs[jobid]['qname']] += Jobs[jobid]['qtime']

    UserCount[Jobs[jobid]['user']]['total'] += 1
    UserRTime[Jobs[jobid]['user']]['total'] += Jobs[jobid]['rtime']
    UserQTime[Jobs[jobid]['user']]['total'] += Jobs[jobid]['qtime']

fig, axs = plt.subplots(2, 2, gridspec_kw={'width_ratios': [1, 1]})
df = pd.DataFrame(UserRTime)
df.sum(axis=1).plot.barh(stacked=False, title='Run Time', ax=axs[0,0], cmap=cmap)
idle=pd.Series({"idle":(dayback*24*108)-df.sum(axis=1)[['total']][0]})
df2=df.sum(axis=1)[['md','fep']].append(idle)
explode = (0.1, 0.1, 0)
df2.plot.pie(title='Run Time', autopct='%1.0f%%', explode=explode, shadow=True, ax=axs[0,1], cmap=cmap)
ax=axs[0,1].set_ylabel('')

df = pd.DataFrame(UserQTime)
df.sum(axis=1).plot.barh(stacked=False, title='Queue Time', ax=axs[1,0], cmap=cmap)
df.sum(axis=1)[['md','fep']].plot.pie(title='Queue Time', shadow=True, ax=axs[1,1], cmap=cmap)
ax=axs[1,1].set_ylabel('')

fig, axs = plt.subplots(2, 2, gridspec_kw={'width_ratios': [1, 1]})
dft = pd.DataFrame(UserRTime).transpose()
colors = []
for ii in dft.sort_values(by=['md']).index:
    colors.append(Colors[ii])
dft.sort_values(by=['md'])[['md']].plot.barh(stacked=False, title='Run Time', ax=axs[0,0], color=[colors])
colors = []
for ii in dft.sort_values(by=['fep']).index:
    colors.append(Colors[ii])
dft.sort_values(by=['fep'])[['fep']].plot(kind="barh", stacked=False, title='Run Time', ax=axs[0,1], color=[colors])

dft = pd.DataFrame(UserQTime).transpose()
colors = []
for ii in dft.sort_values(by=['md']).index:
    colors.append(Colors[ii])
dft[['md']].sort_values(by=['md']).plot(kind="barh", stacked=False, title='Queue Time', ax=axs[1,0], color=[colors])
colors = []
for ii in dft.sort_values(by=['fep']).index:
    colors.append(Colors[ii])
dft[['fep']].sort_values(by=['fep']).plot(kind="barh", stacked=False, title='Queue Time', ax=axs[1,1], color=[colors])

plt.tight_layout()
plt.show()
# plt.savefig('PLOT_pie.png')
