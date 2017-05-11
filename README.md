# Microbial community timeseries analysis

Analysis of persistent states and state transitions in microbial communities from timeseries and static data.

## Authors
* [John D. Chodera](http://choderalab.org) `<john.chodera@choderalab.org>`
* [Libusha Kelly](http://www.einstein.yu.edu/faculty/13827/libusha-kelly/) `<libusha.kelly@einstein.yu.edu>`

## Manifest
* `data/`
    * `BactNorm.txt` - Polz lab Nahant coastal time series (unpublished!)
    * `bats_orig.txt` - microbial sampling from the BATS ocean site
    * `cholera/` - patient microbiome recovery after cholera infection
    * `cogs/` - static samples from 28 ocean sites
    * `david/` - year-long time series of 2 healthy individuals from David 2014
    * `caporaso/` - time series of various body sites and people from Caporaso 2011
* `figures/`
* `notebooks/` - exploratory / explanatory R or Jupyter notebooks
    * `hierarchical-states/` - coarse graining time points with metrics
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

### Data
* David, Lawrence A., Arne C. Materna, Jonathan Friedman, Maria I. Campos-Baptista, Matthew C. Blackburn, Allison Perrotta, Susan E. Erdman, and Eric J. Alm. “Host Lifestyle Affects Human Microbiota on Daily Timescales.” Genome Biology 15 (2014): R89. doi:10.1186/gb-2014-15-7-r89.
* Caporaso, J. Gregory, Christian L. Lauber, Elizabeth K. Costello, Donna Berg-Lyons, Antonio Gonzalez, Jesse Stombaugh, Dan Knights, et al. “Moving Pictures of the Human Microbiome.” Genome Biol 12, no. 5 (2011): R50.
* Hsiao, Ansel, A. M. Shamsir Ahmed, Sathish Subramanian, Nicholas W. Griffin, Lisa L. Drewry, William A. Petri, Rashidul Haque, Tahmeed Ahmed, and Jeffrey I. Gordon. “Members of the Human Gut Microbiota Involved in Recovery from Vibrio Cholerae Infection.” Nature 515, no. 7527 (November 20, 2014): 423–26. doi:10.1038/nature13738.

### Methods
#### tICA
*  Schwantes, Christian R., and Vijay S. Pande. [Improvements in Markov State Model Construction Reveal Many Non-Native Interactions in the Folding of NTL9](http://dx.doi.org/10.1021/ct300878a) J. Chem Theory Comput. 9.4 (2013): 2000-2009.
*  Perez-Hernandez, Guillermo, et al. [Identification of slow molecular order parameters for Markov model construction](http://dx.doi.org/10.1063/1.4811489) J Chem. Phys (2013): 015102.
*  Naritomi, Yusuke, and Sotaro Fuchigami. [Slow dynamics in protein fluctuations revealed by time-structure based independent component analysis: The case of domain motions](http://dx.doi.org/10.1063/1.3554380) J. Chem. Phys. 134.6 (2011): 065101.

#### MSMs
* Prinz, J.-H., et al. [Markov models of molecular kinetics: Generation and validation](http://dx.doi.org/10.1063/1.3565032>) J. Chem. Phys. 134.17 (2011): 174105.
* Pande, V. S., K. A. Beauchamp, and G. R. Bowman. [Everything you wanted to know about Markov State Models but were afraid to ask](http://dx.doi.org/10.1016/j.ymeth.2010.06.002) Methods 52.1 (2010): 99-105.

#### information theory and Jensen-Shannon divergence
* Jianhua Lin. Divergence measures based on the Shannon entropy. IEEE Transactions on Information theory, 37(1):145–151, 1991.
* Omry Koren, Dan Knights, Antonio Gonzalez, Levi Waldron, Nicola Segata, Rob Knight, Curtis Huttenhower, and Ruth E. Ley. A Guide to Enterotypes across the Human Body: Meta-Analysis of Microbial Community Structures in Human Microbiome Datasets. PLOS Computational Biology, 9(1):e1002863, January 2013.
* David, Lawrence A., Arne C. Materna, Jonathan Friedman, Maria I. Campos-Baptista, Matthew C. Blackburn, Allison Perrotta, Susan E. Erdman, and Eric J. Alm. “Host Lifestyle Affects Human Microbiota on Daily Timescales.” Genome Biology 15 (2014): R89. doi:10.1186/gb-2014-15-7-r89.

