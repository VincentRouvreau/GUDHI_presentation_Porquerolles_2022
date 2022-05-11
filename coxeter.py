# # The equations for the 8-knot
#
# ## 1. The knot
#
# The 3-sphere in $\mathbb{R}⁴$ : $S_3 = x² + y² + z² + t² -1$
#
# The 1-sphere in $\mathbb{R}²$ : $S_1 = cos(\alpha)² + sin(\alpha)² -1$
#
# The knot (1) is defined by:
# $$
# \rho = x² + y² + (z² - t²)² + (2zt)²
# $$
# $$
# f_1 = (z² - t²)\rho + x(8x² -2\rho)
# $$
# $$
# f_2 = 2\sqrt{2}xzt + y(8x² -2\rho)
# $$
#
# Then the knot in $\mathbb{R}⁴$ is:
# $$
# S_3 = f_1 = f_2 = 0
# $$
#
# A point on the knot (a valid seed) is the point $X_0 = (x, y, z, t) = (0, 0, 1/\sqrt{2}, 1/\sqrt{2})$
#
# The sources are available [here](coxeter/knot.cpp).

# ### Visualization
# Requires [medit](https://www.ljll.math.upmc.fr/frey/software.html)

# ! medit coxeter/knot.mesh

# Or on this [YouTube link](https://www.youtube.com/shorts/6GdZm9yhAjw)

# ## 2. The fibration of the complement of the knot
#
# In $\mathbb{R}⁴ * \mathbb{R}²$ , the fibration (2) (desingularized along the knot) is defined by:
# $$
# S_1 = S_3 = 0
# $$
# $$
# cos(\alpha)f_1 + sin(\alpha)f_2 = 0
# $$
# $$
# sin(\alpha)f_1 - cos(\alpha)f_2 > 0
# $$
#
# When one fixes $\alpha$ so that $S_1 = 0$ in (2), the output is a Seifert surface whose boundary is the 8-knot. The strict inequality can be replaced by $\geq$ 0.
# We can compute a PL approximation of such a surface by applying the algorithm in $\mathbb{R}⁴$. The output is a 2-dimensional polyhedral surface embedded in $\mathbb{R}⁴$ (in fact in the sphere $S_3 = 0$ of $\mathbb{R}⁴$).
# Note that since any such surface contains the 8-knot, the seed point used for the knot is a valid seed point for any Seifert surface.
# A PL approximation of the entire fibration of dimension 3 can be obtained by running the algorithm in $\mathbb{R}⁶$ using (2).
# We can again replace the strict equality by a non strict one so that the same seed point $X_0 = (x, y, z, t) = (0, 0, 1/\sqrt{2}, 1/\sqrt{2})$ is still valid.
#
# The sources are available  [here](coxeter/seifert_surface.cpp).

# ### Visualization
#
# #### When $\alpha = 0$:

# ! medit coxeter/seifert_surface.mesh

# Or on this [YouTube link](https://www.youtube.com/watch?v=YJ1706hJzmA)
#
# #### For every other $\alpha$:
#
# Please refer to [this webapp](https://mybinder.org/v2/gh/gudhi/coxeter_triangulation_over_fibration_visualization/main?urlpath=%2Fvoila%2Frender%2FSeifert_surface_visualization.ipynb) (quite long to launch as it installs required visualization packages).
