#!/usr/bin/env python

"""
Build a Markov state model of microbial community populations.

@author: John D. Chodera <john.chodera@choderalab.org>

"""

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn

import numpy as np
import pandas as pd

filename = 'BATS_ECOTYPE.txt'

# Load data
feature_names = ['Abundance_e9312', 'Abundance_eMED4', 'Abundance_NATL', 'Abundancee_SS120', 'Abundancee_9313']
nfeatures = len(feature_names)
dtypes = { feature : np.float32 for feature in feature_names }
df = pd.read_table(filename, delim_whitespace=True, parse_dates=['date'], dtype=dtypes, na_values=['nd'])

# Extract depths
depths = np.array(list(set(df['depth'])), np.int32)
depths = np.sort(depths)
ndepths = len(depths)
print('ndepths = %d' % ndepths)

#
# Build MSMs
#

import msmbuilder.cluster
import msmbuilder.hmm
import sklearn.cluster
from msmbuilder.decomposition import tICA

n_components = 5 # number of dimensions to retain
lag_time = 1 # lag time (in sampling intervals)

#
# Perform tICA clustering on each depth separately
#

# Get number of samples
depth = 1
features = np.log(np.array(df[df['depth'] == depth][feature_names])) # log abundances
nsamples = features.shape[0]
print('nsamples = %d' % nsamples)


# Get colormap
from matplotlib.pyplot import cm
cmap = cm.rainbow(np.linspace(0,1,nsamples))

# Make plots
pdf_filename = 'BATS_ECOTYPE-tICA.pdf'
with PdfPages(pdf_filename) as pdf:
    for (depth, color) in zip(depths, cmap):
        plt.figure(figsize=(6,6))
        plt.gca().set_color_cycle(cmap)
        try:
            print('depth = %d' % depth)
            features = np.log(np.array(df[df['depth'] == depth][feature_names])) # log abundances
            tica = tICA(n_components=n_components, lag_time=lag_time)
            projection = tica.fit_transform([features])
            print('eigenvalues:')
            print tica.eigenvalues_
            projection = projection[0] # select out only trajectory
            #plt.plot(projection[:,0], projection[:,1], '.-', color=color)
            for i in range(nsamples-1):
                plt.gca().plot(projection[i:i+2,0], projection[i:i+2,1])

            plt.title('depth %d' % depth)
            plt.xlabel('tIC 1')
            plt.ylabel('tIC 2')
            plt.xticks([])
            plt.yticks([])
            #plt.legend(depths)
            pdf.savefig()
            plt.close()
        except Exception as e:
            print(e)
