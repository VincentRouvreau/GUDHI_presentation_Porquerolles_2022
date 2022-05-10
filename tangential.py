import ipyvolume as ipv
# https://gudhi.inria.fr/python/latest/datasets_generators.html#points-on-sphere
from gudhi.datasets.generators import points

sph = points.sphere(n_samples = 1500, ambient_dim = 3, radius = 1., sample = 'random')
# Display
fig = ipv.figure()
scatter = ipv.scatter(sph[:, 0], sph[:, 1], sph[:, 2], marker="sphere")
ipv.show()

import numpy as np
# https://gudhi.inria.fr/python/latest/point_cloud.html#gudhi.subsampling.choose_n_farthest_points
from gudhi.subsampling import choose_n_farthest_points

sps_sph = choose_n_farthest_points(points=sph, nb_points=1000)
sps_sph = np.array(sps_sph)

fig = ipv.figure()
scatter = ipv.scatter(sps_sph[:, 0], sps_sph[:, 1], sps_sph[:, 2], marker="sphere")
ipv.show()

# https://gudhi.inria.fr/python/latest/tangential_complex_user.html
from gudhi import TangentialComplex

tc = TangentialComplex(intrisic_dim = 2, points=sps_sph)
tc.compute_tangential_complex()
print(f"Tangential contains {tc.num_vertices()} vertices, {tc.num_simplices()} simplices and {tc.num_inconsistent_simplices()} inconsistencies")

st = tc.create_simplex_tree()
print(f"Simplex tree contains {st.num_vertices()} vertices, {st.num_simplices()} simplices and of dimension {st.dimension()}")

# +
tri=[]
for simplex in st.get_skeleton(2):
    if len(simplex[0]) == 3:
        tri.append(simplex[0])

print(f"Simplex tree contains {len(tri)} triangles")
# -

ipv.figure()
mesh = ipv.plot_trisurf(sps_sph[:, 0], sps_sph[:, 1], sps_sph[:, 2], triangles=tri, color='blue')
ipv.scatter(sps_sph[:, 0], sps_sph[:, 1], sps_sph[:, 2], marker='sphere')
ipv.show()
