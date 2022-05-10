# # Rips complex

from gudhi import read_points_from_off_file
# from https://github.com/Ripser/ripser/blob/master/examples/o3_4096.txt
pts = read_points_from_off_file(off_file='data/o3_4096.off')

import profiling_tools as pt

# +
# https://gudhi.inria.fr/python/latest/rips_complex_user.html
from gudhi import RipsComplex

ram = pt.Ram()
chr = pt.Chrono()
rips = RipsComplex(points = pts, max_edge_length=1.4)
print(f"Rips construction tooks {chr.end_timer():.3f} seconds - {ram.get()} Mo RAM")
# -

chr = pt.Chrono()
st = rips.create_simplex_tree(max_dimension=1)
print(f"Number of simplices = {st.num_simplices()} - {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

#  ^
# /!\ This will crash for dimension >= 5 !!
# ---
for dim in range(2, 5):
    chr = pt.Chrono()
    st = rips.create_simplex_tree(max_dimension=dim)
    print(f"Number of simplices = {st.num_simplices()} for dim {dim} - {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

# The number of simplices grows exponentially with the dimension. The ratio *memory consumption / size of complex* is very stable for the simplex tree.

# # Edge collapse
#
# Collapsing edges removes edges that will not modify the persistence results.
#
# ## The importance of working on edges

# Some cleanup
del st
ram = pt.Ram()

chr = pt.Chrono()
st = rips.create_simplex_tree(max_dimension=1)
print(f"Number of simplices = {st.num_simplices()} - {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

chr = pt.Chrono()
st.collapse_edges()
print(f"Number of simplices = {st.num_simplices()} - {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

chr = pt.Chrono()
st.expansion(9)
print(f"Number of simplices = {st.num_simplices()} - {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

# Edge collapse saves memory use for the data structure and computation time. It will also makes persistence faster to be calculted.

chr = pt.Chrono()
diag = st.persistence(persistence_dim_max=True)
print(f"compute persistence tooks {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

# ## Faster - Higher - Stronger
#
# Edge collapse can be launched several times

# Some cleanup
del st
ram = pt.Ram()


# +
def collapse_until(stree):
    nbs = stree.num_simplices() + 1
    while nbs > stree.num_simplices():
        nbs = stree.num_simplices()
        chr = pt.Chrono()
        stree.collapse_edges()
        print(f"Number of simplices = {st.num_simplices()} - {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

st = rips.create_simplex_tree(max_dimension=1)
collapse_until(st)
# -

chr = pt.Chrono()
st.expansion(9)
print(f"Number of simplices = {st.num_simplices()} - {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

chr = pt.Chrono()
diag = st.persistence(persistence_dim_max=True)
print(f"compute persistence tooks {chr.end_timer():.3f} seconds - {ram.get():.1f} Mo RAM")

# +
import matplotlib.pyplot as plt
from gudhi import plot_persistence_diagram

plot_persistence_diagram(diag, legend=True, max_intervals=len(diag))
plt.show()
