import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import MinMaxScaler

import gudhi as gd

D1 = np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.], [0., np.inf], [5., np.inf]])
ax = gd.plot_persistence_diagram(D1)
ax.set_title("D1")
plt.show()

# ## Preprocessing

from gudhi.representations import DiagramSelector, DiagramScaler, Clamping
# https://gudhi.inria.fr/python/latest/representations.html#module-gudhi.representations.preprocessing

proc1 = DiagramSelector(use=True, point_type="finite")
proc2 = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())])
proc3 = DiagramScaler(use=True, scalers=[([1], Clamping(maximum=.9))])
D1 = proc3(proc2(proc1(D1)))

ax = gd.plot_persistence_diagram(D1)
ax.set_title("Finite values of D1 rescaled in [0,1] and less than 0.9")
plt.show()

# ## Vector methods

from gudhi.representations import Landscape, Silhouette, PersistenceImage
# https://gudhi.inria.fr/python/latest/representations.html#module-gudhi.representations.vector_methods

LS = Landscape(resolution=1000)
L = LS(D1)
plt.plot(L[:1000])
plt.plot(L[1000:2000])
plt.plot(L[2000:3000])
plt.title("Landscape")
plt.show()

SH = Silhouette(resolution=1000, weight=lambda x: np.power(x[1]-x[0],2))
plt.plot(SH(D1))
plt.title("Silhouette")
plt.show()

PI = PersistenceImage(bandwidth=.1, weight=lambda x: x[1], im_range=[0,1,0,1], resolution=[100,100])
plt.imshow(np.flip(np.reshape(PI(D1), [100,100]), 0))
plt.title("Persistence Image")
plt.show()


# ## Metrics

D2 = np.array([[1.,5.],[3.,6.],[2.,7.]])
D2 = proc3(proc2(proc1(D2)))

D = [(0, pers) for pers in D1]
D.extend([(1, pers) for pers in D2])

ax = gd.plot_persistence_diagram(D)
ax.set_title("D1 (in red) and D2 (in blue)")
plt.show()

from gudhi.representations import SlicedWassersteinDistance, SlicedWassersteinKernel, WassersteinDistance, BottleneckDistance

swd = SlicedWassersteinDistance(num_directions=100)
print(f"SW distance is {swd(D1, D2)}")

swk = SlicedWassersteinKernel(num_directions=100, bandwidth=1.)
print(f"SW kernel is {swk(D1, D2)}")

wdp = WassersteinDistance(order=2, internal_p=2, mode="pot")
print(f"Wasserstein distance (POT) is {wdp(D1, D2)}")

wdh = WassersteinDistance(order=2, internal_p=2, mode="hera", delta=0.0001)
print(f"Wasserstein distance (hera) is {wdh(D1, D2)}")

bd = BottleneckDistance(epsilon=.001)
print(f"Bottleneck distance is {bd(D1, D2)}")

# ### [Perslay](https://github.com/MathieuCarriere/tda-tutorials/blob/perslay/Tuto-GUDHI-perslay-visu.ipynb)

# ## Scikit-learn like interfaces - *Work in progress*

from sklearn.svm import SVC
from sklearn.pipeline import Pipeline

from gudhi.sklearn.rips_persistence import RipsPersistence
from gudhi.representations import Landscape, DiagramSelector

pipe = Pipeline(
    [
        ("rips_pers", RipsPersistence(max_rips_dimension=2, max_persistence_dimension=2, only_this_dim=1, n_jobs=-1)),
        ("finite_diags", DiagramSelector(use=True, point_type="finite")),
        ("persim", PersistenceImage(bandwidth=50, weight=lambda x: x[1] ** 2, im_range=[0, 256, 0, 256], resolution=[20, 20])),
        ("SVC", SVC()),
    ]
)

# Learn from the train subset
tda_ml_pipe.fit(X_train, y_train)
# Predict from the test subset
predicted = tda_ml_pipe.predict(X_test)
