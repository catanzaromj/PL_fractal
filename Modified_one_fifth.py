"""
This script computes the persistence landscape of the modified Cantor middle fifth set.

authors: Lee Przybylski, Michael J. Catanzaro
"""

import numpy as np
from ripser import ripser
from persim.landscapes import PersLandscapeApprox
from persim.landscapes import plot_landscape
import matplotlib.pyplot as plt
from persim import plot_diagrams

#%% Initialize

Att_prox = np.array([0, 1])


def psi0(x):
    return x / 5


def psi1(x):
    return (x + 1) / 5


def psi4(x):
    return (x + 4) / 5


#%% Compute an approximation of the IFS.
# The scale parameter controls how many iterations of the IFS to compute. The
# larger the scale, the better the approximation, but the more computational
# time.

scale = 5

for _ in range(scale):
    out0 = psi0(Att_prox)
    out1 = psi1(Att_prox)
    out4 = psi4(Att_prox)
    Att_prox = np.concatenate((out0, out1, out4))

# print(np.sort(Att_prox))

l = len(Att_prox)
h = np.zeros(l)

#%% Plot the approximation in the interval
hminus = np.zeros(len(out0))
plt.figure(figsize=(12, 1))
plt.scatter(Att_prox, h, linewidth=0.1)
plt.scatter(out0, hminus, color="red", linewidth=0.1)
plt.scatter(out1, hminus, color="blue", linewidth=0.1)
plt.scatter(out4, hminus, color="green", linewidth=0.1)
plt.scatter(0.2, 0, color="purple", linewidth=0.1)
plt.text(0.5, -0.01, "$S_n$", size=15)
plt.yticks([])

#%% Plot the point in the scale:
out0 = psi0(Att_prox)
out1 = psi1(Att_prox)
out4 = psi4(Att_prox)
Snplus1 = np.concatenate((out0, out1, out4))
hplus = np.zeros(len(Snplus1))
fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 1), sharex=True)
ax1.scatter(Att_prox, h, linewidth=0.1)
# ax1.title = "$S_n$"
ax2.scatter(out0, h, linewidth=0.1, color="red")
ax2.scatter(out1, h, linewidth=0.1, color="blue")
ax2.scatter(out4, h, linewidth=0.1, color="green")
ax2.scatter(0.2, 0, linewidth=0.1, color="purple")
ax2.text(0.5 / 5, -0.01, "$K_1$", size=15)
ax2.text(0.31, -0.01, "$K_2$", size=15)
ax2.text(0.91, -0.01, "$K_3$", size=15)
ax1.text(0.5, -0.01, "$S_n$", size=15)
ax2.text(0.5, -0.01, "$S_{n+1}$", size=15)
ax1.set_yticks([])
ax2.set_yticks([])

#%% Compute the persistent homology for the set
Att_prox = np.reshape(Att_prox, [l, 1])
Dgmprox = ripser(Att_prox, maxdim=0)["dgms"]
plot_diagrams(Dgmprox, show=True, lifetime=True)
#%% Compute the persistence landscape

Landprox = PersLandscapeApprox(dgms=Dgmprox, hom_deg=0)
# Landprox.compute_landscape()
ttl = "PL for $H_0\mathrm{Cech}(\mathcal{C}_{1/5})$"
plot_landscape(
    landscape=Landprox,
    title=ttl,
    labels=["", "", ""],
    colormap="default",
)
plt.xticks([])
plt.yticks([])
#%% Plot a simple version of the landscape, this is faster
plot_landscape(landscape=Landprox, depth_range=range(1, 80), num_steps=5000)
