# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 21:55:30 2022

@author: Lee Przybylski
"""
import numpy as np
from ripser import ripser
from persim.landscapes import PersLandscapeApprox
from persim.landscapes import plot_landscape
import matplotlib.pyplot as plt
from persim import plot_diagrams

#This IFS has fixed points that make up vertices of a regular pentagon

#Fixed points:
pi = np.pi
    
x1 = np.array([0,0])
x2 = np.array([1,0])
x3 = np.array([np.cos(3*pi/5), np.sin(2*pi/5)])
x4 = np.array([0.5, np.sin(2*pi/5) + np.sin(pi/5)])
x5 = np.array([1+ np.cos(2*pi/5), np.sin(2*pi/5)])
x6 = np.array([0.5, 0.5*(np.sin(2*pi/5) + np.sin(pi/5))])

fixed = np.array([x1,x2,x3,x4,x5, x6])

#%% Plot the fixed points
plt.figure(figsize = (10,10))
plt.scatter(fixed[:, 0], fixed[:, 1], linewidth=0.1, Color="blue", zorder=2)
plt.title("$S_0$", size = 20)

#%% Define the IFS
m = 1/6

def phi1(x):
    return m*x 


def phi2(x):
    return m*x + (1-m)*x2


def phi3(x):
    return m*x + (1-m)*x3

def phi4(x):
    return m*x + (1-m)*x4

def phi5(x):
    return m*x + (1-m)*x5

def phi6(x):
    return m*x + (1-m)*x6
#%%
SArray = fixed

#Choose a scale to compute to
scale = 2


for _ in range(scale):
    out1 = phi1(SArray)
    out2 = phi2(SArray)
    out3 = phi3(SArray)
    out4 = phi4(SArray)
    out5 = phi5(SArray)
    out6 = phi6(SArray)
    SArray = np.concatenate((out1, out2, out3, out4, out5, out6), axis=0)

#%% Plot at scale
plt.figure(figsize=(10, 10))
# Color the shape
x = SArray[:, 0]
y = SArray[:, 1]
plt.scatter(SArray[:, 0], SArray[:, 1], linewidth=0.1, Color="blue", zorder=2)
plt.title("$S_2$", size=20)
#%% Compute the persistence
Dgmprox = ripser(SArray, maxdim=1)["dgms"]
plot_diagrams(Dgmprox, show=True, lifetime=False)
print(Dgmprox)
#%% Compute the landscape
Landprox = PersLandscapeApprox(dgms=Dgmprox, hom_deg=0)
ttl = f"PL of Cantor triangle at scale {scale}"
plot_landscape(landscape=Landprox, labels=["", "", ""], title=ttl)
