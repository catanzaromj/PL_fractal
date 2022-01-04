"""
This script computes the persistence landscape of the  Cantor triangle.

authors: Lee Przybylski, Michael J. Catanzaro
"""

import numpy as np
from ripser import ripser
from persim.landscapes import PersLandscapeApprox
from persim.landscapes import plot_landscape
import matplotlib.pyplot as plt
from persim import plot_diagrams


#%% Initialize the Cantor triangle example
def phi0(x):
    return x / 3


def phi20(x):
    return (x + np.array([2, 0])) / 3


def phi02(x):
    return (x + np.array([0, 2])) / 3


SArray = np.array([[0, 0], [1, 0], [0, 1]])
#%% Chose a scale to compute the ifs up to
# Keep the scale less than 30
scale = 4


for _ in range(scale):
    out0 = phi0(SArray)
    out20 = phi20(SArray)
    out02 = phi02(SArray)
    SArray = np.concatenate((out0, out20, out02), axis=0)

#%% Plot at scale
plt.figure(figsize=(5, 5))
# Color the shape
x = SArray[:, 0]
y = SArray[:, 1]
plt.scatter(SArray[:, 0], SArray[:, 1], linewidth=0.1, Color="blue", zorder=2)
plt.title("$S_3$", size=20)
#%% Compute the persistence
Dgmprox = ripser(SArray, maxdim=1)["dgms"]
plot_diagrams(Dgmprox, show=True, lifetime=False)
print(Dgmprox)
#%% Compute the landscape
Landprox = PersLandscapeApprox(dgms=Dgmprox, hom_deg=0)
ttl = f"PL of Cantor triangle at scale {scale}"
plot_landscape(landscape=Landprox, labels=["", "", ""], title=ttl)
