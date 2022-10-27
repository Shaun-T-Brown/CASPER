# CASPER (Concentration And Shape Parameter Estimation Routine)

Casper is a python package aimed at predicting the concentration and shape parameter of dark matter haloes as a function of mass and redshift for a specified cosmology.

### Requirements

The module requires the following:

- python3.6 or higher

### Installation

The easiest way to install Casper is using pip:

```
pip install py-casper [--user]
```

The --user flag may be required if you do not have root privileges. Alternatively for a more involved 'installation' that is also editable you can simply  clone the github repository:

```
git clone https://github.com/Shaun-T-Brown/CASPER.git
```

and add the folder (specifically src) to your python path. This method will also require scipy and numpy to be installed. Using pip will automatically install all dependencies.

To use the main functionality of the package (i.e. to predict c and alpha) you will need to be able to generate the linear power spectra for a given cosmology. In principle this can be done using any reliable method. However, we recommend installing and using [CAMB](https://camb.readthedocs.io/en/latest/).



### Usage

The best way to demonstrate how Casper can be used is with a few examples. An interactive jupyter notebook with all necessary modules available can be found on [Binder](https://mybinder.org/v2/gh/Shaun-T-Brown/CASPER-example.git/HEAD?filepath=.%2Fexample_script.ipynb). The Binder notebook can sometimes take a long time to load, a static version of the same notebook can be found at the following github [repository](https://github.com/Shaun-T-Brown/CASPER-example.git), specifically in the file *example_script.ipynb*.

The functionality of the package is relatively straightforward and is essentially a collection of useful functions centered around the density profiles of dark matter haloes. Basic documentation can be generated using *pydoc*.


```
import casper

help(casper)
```


### Acknowledging the code

If the results of this code, particularly the predictions for halo concentration and the shape parameter, are used in any published work then please acknowledge and cite the original [paper](https://arxiv.org/abs/2110.01632) appropriately.

The following bibtex entry may be used:

```
@ARTICLE{2022MNRAS.509.5685B,
       author = {{Brown}, Shaun T. and {McCarthy}, Ian G. and {Stafford}, Sam G. and {Font}, Andreea S.},
        title = "{Towards a universal model for the density profiles of dark matter haloes}",
      journal = {\mnras},
     keywords = {methods: numerical, cosmology: theory, dark matter, Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Astrophysics of Galaxies},
         year = 2022,
        month = feb,
       volume = {509},
       number = {4},
        pages = {5685-5701},
          doi = {10.1093/mnras/stab3394},
archivePrefix = {arXiv},
       eprint = {2110.01632},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.5685B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

For any questions and enquires please contact me via email at *Shaun.T.Brown@durham.ac.uk*


