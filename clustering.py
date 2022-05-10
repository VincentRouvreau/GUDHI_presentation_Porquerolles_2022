# # Clustering
# In this notebook, we will go through an easy example to present ToMaTo, a persistance based clustering tool.

# +
import matplotlib.pyplot as plt

cmap = plt.cm.Spectral;
fig, ax = plt.subplots();

import random as rd
import numpy as np


# Simple function to get random values for x uniformly but within intervals (0,a) U (b, 1)
def x_var(x):
    if x > 0.5:
        return rd.uniform(0.6, 1)
    else:
        return rd.uniform(0, 0.4)

    
p1 = np.zeros((200,2))
for i in range(200):
    p1[i,0] = x_var(rd.uniform(0,1))
    p1[i,1] = rd.uniform(0,1)
    
ax.cla()
ax.scatter(*zip(*p1))
# -

# ## Introduction
#
# This code is an implemention of the ToMATo algorithm, a clustering method based on the idea of topological persistance. In short, the algorithm needs a density estimation (so to each point $x$ we associate a value $\hat{f}(x)$) and a neighborhood graph.
#
# ## Density estimation
#
# ### KDE (from scikit-learn)

# +
from sklearn.neighbors import KernelDensity

kde = KernelDensity(bandwidth=0.1).fit(p1)
ax = plt.axes(projection='3d')
z_kde = kde.score_samples(p1)
ax.scatter(p1[:, 0], p1[:, 1], z_kde, c=z_kde)
# -

# https://gudhi.inria.fr/python/latest/clustering.html
from gudhi.clustering.tomato import Tomato

clust = Tomato(
        density_type = 'KDE',
        r = 0.1, # aka. bandwidth
        # ...
    )

# ### DTM

# https://gudhi.inria.fr/python/latest/point_cloud.html#gudhi.point_cloud.dtm.DTMDensity
from gudhi.point_cloud.dtm import DTMDensity

z_dtm = DTMDensity(k=30, normalize=True).fit_transform(p1)
ax = plt.axes(projection='3d')
ax.scatter(p1[:, 0], p1[:, 1], z_dtm, c=z_dtm)

clust = Tomato(
        density_type = 'DTM',
        k = 30,
        # ...
    )

# ## Neighborhood graph
#
# ### K-nearest neighbors with $k = 4$

# https://gudhi.inria.fr/python/latest/point_cloud.html#module-gudhi.point_cloud.knn
from gudhi.point_cloud.knn import KNearestNeighbors

# +
X = np.array(p1)
nbrs = KNearestNeighbors(k=4)
indices = nbrs.fit_transform(X)
plt.plot(X[:,0], X[:,1], 'o')
for i in indices:
    Y = np.zeros((2,2))
    for j in range(len(i)):
        Y[0][0]= X[int(i[0])][0]
        Y[1][0]= X[int(i[0])][1]
        Y[0][1]= X[int(i[j])][0]
        Y[1][1]= X[int(i[j])][1]
        plt.plot(Y[0], Y[1], 'ro-')

plt.show()
# -

clust = Tomato(
        graph_type = 'knn',
        k = 4,
        # ...
    )

# ### Radius graph with $r = 0.1$

# +
nbrs = KNearestNeighbors(k=30, return_distance= True)
indices, distances = nbrs.fit_transform(X)
plt.plot(X[:,0], X[:,1], 'o')
for i in indices:
    Y = np.zeros((2,2))
    for j in range(len(i)):
        if distances[int(i[0]), j] < 0.1:
            Y[0][0]= X[int(i[0])][0]
            Y[1][0]= X[int(i[0])][1]
            Y[0][1]= X[int(i[j])][0]
            Y[1][1]= X[int(i[j])][1]
            plt.plot(Y[0], Y[1], 'ro-')

plt.show()
# -

clust = Tomato(
        graph_type = 'radius',
        r = 0.1,
        # ...
    )

# ## ToMaTo
#
# Tomato class has been designed to be scikit-learn like kind of interface.

# +
ex1 = Tomato(
        graph_type="radius",
        density_type="KDE",
        r=0.1,
    )

labels = ex1.fit_predict(p1)
print(labels)

print("\nThere are " + str(ex1.n_clusters_) + " initial clusters")
ex1.plot_diagram()
# -

ex1.n_clusters_ = 2
print(ex1.n_clusters_)
print(ex1.labels_)

fig, ax = plt.subplots()
ax.cla()
ax.scatter(*zip(*p1), c=ex1.labels_)

# If the connected components from the neighborhood graph is more than 2, consider increasing the radius to $r = 0.15$ to construct the graph.
