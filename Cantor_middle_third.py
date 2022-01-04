#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script computes the persistence landscape of the Cantor middle third set.

authors: Lee Przybylski, Michael J. Catanzaro
"""

import numpy as np
from ripser import ripser
from persim.landscapes import PersLandscapeApprox
from persim.landscapes import plot_landscape
import matplotlib.pyplot as plt
from persim import plot_diagrams

#%% Initialize the middle third Cantor set
#
#
# In this script, we construct the middle third Cantor set up to a specified
# scale. We proceed using the IFS approach to affine fractals.

C3prox = np.array([0, 1])


def phi0(x):
    return x / 3


def phi2(x):
    return (x + 2) / 3


#%% Compute an approximation of the IFS.
# The scale parameter controls how many iterations of the IFS to compute. The
# larger the scale, the better the approximation, but the more computational
# time.

scale = 8

for _ in range(scale):
    out0 = phi0(C3prox)
    out2 = phi2(C3prox)
    C3prox = np.concatenate((out0, out2))

# print(np.sort(C3prox))

l = int(np.exp2(scale + 1))
h = np.zeros(l)

#%% Plot the approximation in the interval

plt.figure(figsize=(12, 1))
plt.scatter(C3prox, h, linewidth=0.01)
plt.text(0.5, -0.01, "$S_n$", size=15)
plt.yticks([])

out0 = phi0(C3prox)
out2 = phi2(C3prox)
Snplus1 = np.concatenate((out0, out2))
hplus = np.zeros(2 * l)
fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 1), sharex=True)
ax1.scatter(C3prox, h, linewidth=0.1)
# ax1.title = "$S_n$"
ax2.scatter(out0, h, linewidth=0.1, color="green")
ax2.scatter(out2, h, linewidth=0.1, color="red")
# ax2.scatter(Snplus1, hplus, linewidth = .1)
ax2.text(0.16, -0.01, "$L_n$", size=15)
ax2.text(0.82, -0.01, "$R_n$", size=15)
ax1.text(0.5, -0.01, "$S_n$", size=15)
# ax2.text(0.5, -.01, "$S_{n+1}$", size= 15)
ax1.set_yticks([])
ax2.set_yticks([])

#%% Compute the zeroth persistent homology of the approximaation

C3prox = np.reshape(C3prox, [l, 1])
Dgmprox = ripser(C3prox, maxdim=0)["dgms"]
plot_diagrams(Dgmprox, show=True, lifetime=True)

#%% Compute the persistence landscape

Landprox = PersLandscapeApprox(dgms=Dgmprox, hom_deg=0)
ttl = "PL for $H_0\mathrm{Cech}(\mathcal{C})$"
plot_landscape(
    landscape=Landprox,
    title=ttl,
    labels=["", "", ""],
    colormap="default",
)
plt.title("PL for $H_0\mathrm{Cech}(\mathcal{C})$", size=20)
plt.xticks([])
plt.yticks([])
