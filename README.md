# Microbial community timeseries analysis

Analysis of persistent states and state transitions in microbial communities from timeseries and static data.

## Authors
* [John D. Chodera](http://choderalab.org) `<john.chodera@choderalab.org>`
* [Libusha Kelly](http://www.einstein.yu.edu/faculty/13827/libusha-kelly/) `<libusha.kelly@einstein.yu.edu>`

## Manifest
* `data/`
    * `BactNorm.txt`
    * `bats_orig.txt` - microbial sampling from the BATS ocean site
    * `cholera/` - patient microbiome recovery after cholera infection
    * `cogs/` - static samples from 28 ocean sites
* `figures/`
* `notebooks/` - exploratory / explanatory R or Jupyter notebooks
* `references/` - PDFs
* `scripts/` - code to analyze `data/` and generate `figures/`

## Installing dependencies

Install [miniconda](http://conda.pydata.org/miniconda.html), then dependencies:
On `osx` using `bash`:
```bash
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
bash ./Miniconda2-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda2
export PATH="$HOME/miniconda2/bin:$PATH"
conda config --add channels omnia
conda install --yes numpy pandas pyemma seaborn scikit-learn msmbuilder
```


## References

### tICA
*  Schwantes, Christian R., and Vijay S. Pande. [Improvements in Markov State Model Construction Reveal Many Non-Native Interactions in the Folding of NTL9](http://dx.doi.org/10.1021/ct300878a) J. Chem Theory Comput. 9.4 (2013): 2000-2009.
*  Perez-Hernandez, Guillermo, et al. [Identification of slow molecular order parameters for Markov model construction](http://dx.doi.org/10.1063/1.4811489) J Chem. Phys (2013): 015102.
*  Naritomi, Yusuke, and Sotaro Fuchigami. [Slow dynamics in protein fluctuations revealed by time-structure based independent component analysis: The case of domain motions](http://dx.doi.org/10.1063/1.3554380) J. Chem. Phys. 134.6 (2011): 065101.

### MSMs
* Prinz, J.-H., et al. [Markov models of molecular kinetics: Generation and validation](http://dx.doi.org/10.1063/1.3565032>) J. Chem. Phys. 134.17 (2011): 174105.
* Pande, V. S., K. A. Beauchamp, and G. R. Bowman. [Everything you wanted to know about Markov State Models but were afraid to ask](http://dx.doi.org/10.1016/j.ymeth.2010.06.002) Methods 52.1 (2010): 99-105.
