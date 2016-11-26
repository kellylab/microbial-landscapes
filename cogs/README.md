# Static microbiome samples at various ocean sites

All Prochlorococcus and viral gene cluster IDs and their length-normalized abundances at 28 ocean sites.

The first column is the name of the gene cluster (`COG0`)

All other columns are sampling sites, they are named as follows:
* `BS1D` = bacterial fraction, station 1, deep
* `PS1S` = phage fraction, station 1, surface

## Manifest

* `BIGRAPA_COGS.csv` - length-normalized abundance ocean sampling data
* `spectral clustering experiments.ipynb` - spectral clustering experiments IPython notebook

## Methods

This example uses [spectral clustering](http://link.springer.com/article/10.1007/s11634-013-0134-6) on a transition matrix constructed from a simple similarity metric.
