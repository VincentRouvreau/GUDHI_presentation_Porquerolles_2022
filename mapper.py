# # Topological Exploratory Data Analysis
# In this notebook, we will see how to perform topological dimension reduction: we will compute simplicial complex approximations of a point cloud. This complex will be a Mapper complex. This complex uses covers of the initial space (such as Voronoi partitions or preimages of filter functions), and use these covers to generate simplicial complexes, by taking the nerve (Mapper).

from gudhi.sklearn import CoverComplex
import numpy as np
import matplotlib.pyplot as plt
import ipyvolume as ipv

X = np.loadtxt('data/human2.txt')

# Display
fig = ipv.figure()
scatter = ipv.scatter(X[::2,1], X[::2,0], X[::2,2], marker="sphere", color="blue")
ipv.show()

# We will use the height function to color the complex nodes.

height = X[:,2]

# Mapper complex with a preimage cover with automatic resolution and hierarchical clustering obtained with automatic threshold.

cover_complex = CoverComplex(
    complex_type='mapper', input_type='point cloud', cover='functional', colors=height[:,np.newaxis], mask=0,
    clustering=None, N=100, beta=0., C=10,
    filters=height[:,np.newaxis], filter_bnds=None, resolutions=None, gains=None,
    input_name="human", cover_name="coord2", color_name="coord2", verbose=True)

# ## Complex computation

_ = cover_complex.fit(X)

# ## Visualization
#
# You can visualize the complex with networkx:

import networkx as nx
G = cover_complex.get_networkx()
plt.figure()
color = [cover_complex.node_info[v]["colors"][0] for v in G.nodes()]
size = [cover_complex.node_info[v]["size"] for v in G.nodes()]
nx.draw(G, pos=nx.kamada_kawai_layout(G), node_color=color, node_size=size)
plt.show()
