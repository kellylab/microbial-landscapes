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
print('depths:')
print(depths)
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

# Only plot some depths
depths = [1, 10, 20, 40, 120, 140, 160, 180]

# Make plots
subplot_index = 0 # index of current subplot
subplot_ny = 2 # number of plots in y direction
subplot_nx = 4 # number of plots in x direction
pdf_filename = 'BATS_ECOTYPE-tICA.pdf'
with PdfPages(pdf_filename) as pdf:
    fig = plt.figure(figsize=(7.5,2.5))
    for (depth, color) in zip(depths, cmap):
        # Select subplot
        subplot_index += 1 # increment subplot index
        ax = plt.subplot(subplot_ny,subplot_nx,subplot_index)
        ax.set_color_cycle(cmap)
        try:
            # Extract features for this depth
            print('depth = %d' % depth)
            features = np.log(np.array(df[df['depth'] == depth][feature_names])) # log abundances

            # Try to compute tICA projection
            # NOTE: This sometimes failes because of missing values, so this is in a try...except clause
            tica = tICA(n_components=n_components, lag_time=lag_time)
            projection = tica.fit_transform([features])

            # Print eigenvalues
            print('eigenvalues:')
            print tica.eigenvalues_
            projection = projection[0] # select out only trajectory

            # Plot with time-varying color
            for i in range(nsamples-1):
                ax.plot(projection[i:i+2,0], projection[i:i+2,1])

            # Title subplot
            plt.title('%d m' % depth, fontsize=9)

            # Fix up axes.
            if (subplot_index-1) % subplot_nx == 0:
                plt.ylabel('tIC 2', fontsize=7)
            if (subplot_index-1) / subplot_nx == 1:
                plt.xlabel('tIC 1', fontsize=7)
            plt.xticks([])
            plt.yticks([])

            if subplot_index % subplot_nx == 0:
                ax2 = ax.twinx()
                if (subplot_index <= subplot_nx):
                    ax2.set_ylabel('shallow', fontsize=7, weight='bold', rotation=90)
                else:
                    ax2.set_ylabel('deep', fontsize=7, weight='bold', rotation=90)
                ax2.set_xticks([])
                ax2.set_yticks([])


        except Exception as e:
            print(e)

    # Write figure
    pdf.savefig()
    plt.close()
