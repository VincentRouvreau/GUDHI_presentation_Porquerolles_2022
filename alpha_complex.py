# # Alpha (and Delaunay) complex
#
# If your data are in an Euclidean space and in a low dimension ($<7$).

import numpy as np
# https://gudhi.inria.fr/python/latest/datasets_generators.html#points-on-sphere
from gudhi.datasets.generators import points


def N_sphere(n_samples, ambient_dim, radius, noise, origin = None):
    sph = points.sphere(n_samples = n_samples, ambient_dim = ambient_dim, radius = radius, sample = 'random')
    sph = sph + noise * np.random.rand(sph.shape[0],sph.shape[1])
    if origin is not None:
        sph = sph + origin * np.ones((sph.shape[0],sph.shape[1]))
    return sph


# sph4 is a 4-sphere
sph4 = N_sphere(n_samples = 500, ambient_dim = 5, radius = 10., noise = 1.)
# sph3 is a 3-sphere
sph3 = N_sphere(n_samples = 1000, ambient_dim = 4, radius = 7., noise = 0.5, origin = [20., 0., 0., 0.])
sph3 = np.insert(sph3, 4, values=0, axis=1)
# concatenate sph3 and sph4
sph = np.concatenate((sph3, sph4), axis=0)

import matplotlib.pyplot as plt
# https://gudhi.inria.fr/python/latest/alpha_complex_user.html
import gudhi as gd

alpha = gd.AlphaComplex(points = sph)
simplex_tree = alpha.create_simplex_tree()
diag = simplex_tree.persistence()
gd.plot_persistence_diagram(diag, legend=True)
plt.show()
