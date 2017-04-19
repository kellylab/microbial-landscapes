# 16S and 18S (relative) normalized counts of microbial taxa (rows, called "OTUs") from ~3 months of daily
# sampling off the coast of Nahant, MA.
# High resolution.

# Files: BactNorm.txt is bacteria, EukNorm.txt is eukaryotes.

# Includes abundance data for large numbers of OTUs..

# Alright, let's take a quick look at this...

import numpy as np
from msmbuilder.decomposition import tICA
from sklearn.decomposition import PCA
import pyemma

def load_norm_data(data_path):
    """This was difficult for pandas to parse, since it's whitespace-delimited but with a
    variable number of apparent columns (since there's whitespace in the lineage strings).

    Pandas throws this error, for example:
    >> df = pd.read_table(data_path, delim_whitespace=True)
    Error: pandas.io.common.CParserError: Error tokenizing data. C error: Expected 92 fields in line 4, saw 93
    """

    # first, read the header
    with open(data_path, "rb") as f:
        header = f.readline().split() # ["OTU", "204", "205", ..., "296", "ConsensusLineage"]

    n_features = len(header) - 1

    # next, read and split all lines of the file after the header
    with open(data_path, "rb") as f:
        lines = [line.split() for line in f.readlines()[1:]]

    # next, take the numerical features from this file and put them in a numpy array
    X = np.array([[float(number) for number in line[1:n_features]] for line in lines])
    print(X.shape)

    # and the OTUs
    OTUs = [line[0] for line in lines]

    # and the consensus lineages
    lineages = [line[n_features:] for line in lines]

    print(header)
    print(len(header))
    print(len(OTUs), OTUs[:10]) # ~30,000 species...
    print(lineages[0])

    return {"header": header, "lineages": lineages, "X": X, "OTUs": OTUs}

# load the datasets
bact_norm = load_norm_data("../data/BactNorm.txt")
euk_norm = load_norm_data("../data/EukNorm.txt")

# They're not sampled at exactly the same time points, but maybe we can just look at
# the overlapping time-points?
print(bact_norm["header"][:-1] == euk_norm["header"][:-1])
print(len(set(bact_norm["header"][:-1]).intersection(set(euk_norm["header"][:-1]))))

# So there are ~90 overlapping points

# Let's see how well we can apply PCA to one matrix or another...

pca = PCA()
pca.fit(bact_norm["X"])
print(np.cumsum(pca.explained_variance_ratio_)[:10])

pca = PCA()
pca.fit(euk_norm["X"])
print(np.cumsum(pca.explained_variance_ratio_)[:10])

# what if we do tICA on a small-ish subset of bacterial species...
tica = tICA()
X = bact_norm["X"]
X_ = X[::10].T
print(X_.shape)
tica.fit([X_])
print("MSMBuilder tICA timescales", tica.timescales_)
print("MSMBuilder tICA eigenvalues", tica.eigenvalues_)
# that's weird: eigenvalues > 1, negative leading timescales..

# does this also happen with pyemma's implementation? this is weird
tica = pyemma.coordinates.tica(X_, lag=1)
print("PyEMMA tICA timescales", tica.timescales)
print("PyEMMA tICA eigenvalues", tica.eigenvalues)

# hmm, I think this is because MSMbuilder uses a shrinkage estimator for the covariance / time-lagged
# covariance matrices...

# let's re-do this with the full set of species?
#X_ = X.T
#tica = pyemma.coordinates.tica(X_, lag=1)

# Okay, cannot do this directly on laptop -- runs out of of memory because EVD of 30k x 30k matrices is hard
# -- can do on cluster or use a transpose trick or use hierarchical tICA...
