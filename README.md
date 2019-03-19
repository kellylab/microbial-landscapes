# Inferring the quasipotential landscapes of microbiomes with topological data analysis

Analysis of persistent states and state transitions in microbial communities from time-series data.

## Authors
* [Libusha Kelly](http://www.einstein.yu.edu/faculty/13827/libusha-kelly/) `<libusha.kelly@einstein.yu.edu>`
* William K. Chang `<william.chang@einstein.yu.edu>`

## Manifest
* `code/` - code to analyze `data/` and generate `figures/`
* `data/`
    * `bats_orig.txt` - microbial sampling from the BATS ocean site
    * `cholera/` - patient microbiome recovery after cholera infection
    * `david/` - year-long time series of 2 healthy individuals from David 2014
    * `hot_orig.txt` - microbial sampling from the HOT ocean site
* `figures/`
* `references/` - PDFs

## Installing dependencies

[Rstudio](http://rstudio.com)

Installing required R packages:

```r
install.packages(c("cowplot", "data.table", "ggraph", "igraph", "philentropy", "tidygraph", tidyverse"))
devtools::install_github("wkc1986/TDAmapper")
```


## References

### Data
- David, Lawrence A., Arne C. Materna, Jonathan Friedman, Maria I. Campos-Baptista, Matthew C. Blackburn, Allison Perrotta, Susan E. Erdman, and Eric J. Alm. “Host Lifestyle Affects Human Microbiota on Daily Timescales.” Genome Biology 15 (2014): R89. doi:10.1186/gb-2014-15-7-r89.
- Hsiao, Ansel, A. M. Shamsir Ahmed, Sathish Subramanian, Nicholas W. Griffin, Lisa L. Drewry, William A. Petri, Rashidul Haque, Tahmeed Ahmed, and Jeffrey I. Gordon. “Members of the Human Gut Microbiota Involved in Recovery from Vibrio Cholerae Infection.” Nature 515, no. 7527 (November 20, 2014): 423–26. doi:10.1038/nature13738.
- Malmstrom, R. R. et al. Temporal dynamics of Prochlorococcus ecotypes in the Atlantic and Pacific oceans. *ISME Journal* 4, 1252–1264 (2010).


### Methods

#### TDA and Mapper
- Rizvi, A. H. et al. Single-cell topological RNA-Seq analysis reveals insights into cellular differentiation and development. *Nature Biotechnology* 35, 551–560. issn: 1087-0156 (June 2017).
- Wasserman, L. Topological Data Analysis. *arXiv*:1609.08227 [stat]. arXiv: 1609.08227. http://arxiv.org/abs/1609.08227 (2018) (Sept. 26, 2016).

#### information theory and Jensen-Shannon divergence
* Jianhua Lin. Divergence measures based on the Shannon entropy. IEEE Transactions on Information theory, 37(1):145–151, 1991.
* Omry Koren, Dan Knights, Antonio Gonzalez, Levi Waldron, Nicola Segata, Rob Knight, Curtis Huttenhower, and Ruth E. Ley. A Guide to Enterotypes across the Human Body: Meta-Analysis of Microbial Community Structures in Human Microbiome Datasets. PLOS Computational Biology, 9(1):e1002863, January 2013.

