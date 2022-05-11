# # DTM Rips complex

import numpy as np
import ipyvolume as ipv
# https://gudhi.inria.fr/python/latest/datasets_generators.html#points-on-sphere
from gudhi.datasets.generators import points

# sph is a 2-sphere
sph = points.sphere(n_samples = 150, ambient_dim = 3, radius = 1., sample = 'random')

outl = (np.random.rand(50,3) * 2.) - 1.

# Display
fig = ipv.figure()
scatter = ipv.scatter(sph[:, 0], sph[:, 1], sph[:, 2], marker="sphere", color="blue")
scatter = ipv.scatter(outl[:, 0], outl[:, 1], outl[:, 2], marker="sphere", color="red")
ipv.show()

# concatenate sphere and outlayers
pts = np.concatenate((sph, outl), axis=0)
print(pts.shape)

# https://gudhi.inria.fr/python/latest/rips_complex_user.html
from gudhi import RipsComplex

rips = RipsComplex(points = pts)
stree = rips.create_simplex_tree(max_dimension=3)
diag = stree.persistence()

# +
import matplotlib.pyplot as plt
from gudhi import plot_persistence_diagram

plot_persistence_diagram(diag, legend=True, max_intervals=len(diag))
plt.show()
# -

# https://gudhi.inria.fr/python/latest/rips_complex_user.html#dtm-rips-complex
from gudhi.dtm_rips_complex import DTMRipsComplex

dtm = DTMRipsComplex(points = pts, k=int(0.05*len(pts)))
stree = dtm.create_simplex_tree(max_dimension=3)
diag = stree.persistence()

# +
import matplotlib.pyplot as plt
from gudhi import plot_persistence_diagram

plot_persistence_diagram(diag, legend=True, max_intervals=len(diag))
plt.show()
