# Markov state models of the temporal dyanmics of microbial ecosystems

Build a Markov state model of microbial community populations.

## Usage

Project the data with [time-structure independent correlation analysis (tICA)](http://msmbuilder.org/3.3.0/tica.html):
```bash
python bats_cluster_tica.py
```

Build a Markov state model of `data/bats_orig.txt`, building separate model for each depth:
```bash
python bats_build_msm.py
```

## Manifest
* `*_cluster_tica.py` - scripts illustrating time-structure independent correlation analysis (tICA)
* `bats_build_msm.py` - simple script illustrating the construction of microbial ecosystem Markov state model
* `cogs-spectral-clustering.py` - loads COGS data, eventually will use [spectral clustering](http://link.springer.com/article/10.1007/s11634-013-0134-6) on a transition matrix constructed from a simple similarity metric?
* `nahant_tica.py` - loads Nahant dataset, attempts to apply tICA directly to OTUs
