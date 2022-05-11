# GUDHI presentation at Porquerolles in 2022

## From a notebook server

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/VincentRouvreau/GUDHI_presentation_Porquerolles_2022/main)

## Installation

It requires some packages to be installed, either:
* with [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)  (recommended, as it is easier to uninstall), then:
```
conda create --name gudhi
conda activate gudhi
conda install -c conda-forge python
pip install -r requirements.txt
```
* with pip (comes with python installation), then `pip install -r requirements.txt`

## Presentation

* Manifold reconstruction:
  * [Tangential complex](tangential.py)
  * [Coxeter triangulation](coxeter.py)

* Simplicial complexes construction:
  * [Alpha complex](alpha_complex.py)
  * [Vietoris Rips complex](rips_complex.py)
  * [DTM filtration](dtm_rips_complex.py)

* Topological descriptors
  * [Representations](persistence_representation.py)

* Clustering:
  * [ToMaTo](clustering.py)

* Cover complex:
  * [Mapper](mapper.py)
