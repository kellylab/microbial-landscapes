#!/usr/bin/env python

"""
Build a Markov state model of microbial community populations from cholera patients.

@author: John D. Chodera <john.chodera@choderalab.org>

"""

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn

import numpy as np
import pandas as pd

filename = 'A_gordon.sorted.txt'

# Load data
df = pd.read_table(filename, delim_whitespace=True)

# Extract sorted sample labels
samples = list(df.index)
nsamples = len(samples)
print('samples = ')
print(samples)
print('nsamples = %d' % nsamples)

# Extract feature names
features = list(df.columns)
nfeatures = len(features)

# Extract feature vectors
patient_timeseries = df.as_matrix()

#
# Analyze
#

import msmbuilder.cluster
import msmbuilder.hmm
import sklearn.cluster
from msmbuilder.decomposition import tICA

n_components = 2 # number of dimensions to retain
lag_time = 1 # lag time (in sampling intervals)

#
# Perform tICA clustering on each depth separately
#

# Get colormap
from matplotlib.pyplot import cm
cmap = cm.rainbow(np.linspace(0,1,nsamples))

# Make plots
pdf_filename = 'cholera.pdf'
patients = ['A']
with PdfPages(pdf_filename) as pdf:
    fig = plt.figure(figsize=(5.5,5.5))
    for patient in patients:
        subplot_index = 1
        subplot_nx = 1
        subplot_ny = 1
        ax = plt.subplot(subplot_ny, subplot_nx, subplot_index)
        ax.set_color_cycle(cmap)
        try:
            # Try to compute tICA projection
            # NOTE: This sometimes failes because of missing values, so this is in a try...except clause
            tica = tICA(n_components=n_components, lag_time=lag_time)
            projection = tica.fit_transform([patient_timeseries])

            # Print eigenvalues
            print('eigenvalues:')
            print tica.eigenvalues_
            projection = projection[0] # select out only trajectory

            # Plot with time-varying color
            for i in range(nsamples-1):
                ax.plot(projection[i:i+2,0], projection[i:i+2,1], 'o-', markersize=3, linewidth=0.5)

            # Title subplot
            plt.title('patient %s' % patient, fontsize=9)

            # Fix up axes.
            plt.ylabel('tIC 2', fontsize=7)
            plt.xlabel('tIC 1', fontsize=7)
            plt.xticks([])
            plt.yticks([])

            if subplot_index % subplot_nx == 0:
                ax2 = ax.twinx()
                ax2.set_xticks([])
                ax2.set_yticks([])


        except Exception as e:
            print(e)

    # Write figure
    pdf.savefig()
    plt.close()
