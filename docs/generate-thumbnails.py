#!/usr/bin/env python
"""
This script generates thumbnails for all the examples 
that deserve one but don't generate one themselves.
"""

import os

import matplotlib.pyplot as plt
from matplotlib.patches import Arc

"""
Thumbnails for Network Generator Examples
=========================================

plot_add_angles.py
------------------
"""

plt.figure()
plt.plot([0, 1], [0, 0], "k-", lw=2)  # Base line
plt.plot([0, 0.5], [0, 1], "k-", lw=2)  # Left line
# add angle marker arc
arc = Arc((0, 0), 0.75, 0.75, angle=0, theta1=0, theta2=62.5, color="k", lw=2)
plt.gca().add_patch(arc)

plt.xlim(-0.2, 1.2)
plt.ylim(-0.2, 1.2)
plt.axis("off")

# Save the sketch as a thumbnail
target_file = os.path.join(
    os.path.dirname(__file__),
    "_static/thumbnails/network_generator/plot_add_angles.png",
)
if not os.path.exists(os.path.dirname(target_file)):
    os.makedirs(os.path.dirname(target_file))
plt.savefig(
    target_file,
    bbox_inches="tight",
    dpi=300,
)
plt.close()

"""
plot_.py
--------------------------
"""
