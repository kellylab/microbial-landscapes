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

def read_patient_data(patient):
    filename = '%s_gordon.taxa_count.sorted.txt' % patient

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

    return [samples, features, patient_timeseries]

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

# Read patient data
patients = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
data = list()
for (subplot_index, patient) in enumerate(patients):
    # Load patient data
    [samples, features, patient_timeseries] = read_patient_data(patient)
    data.append(patient_timeseries)
    nsamples = len(samples)
    nfeatures = len(features)
tica = tICA(n_components=n_components, lag_time=lag_time)
projected_data = tica.fit_transform(data)

# Make plots
pdf_filename = 'cholera.pdf'
subplot_nx = 2
subplot_ny = 4
with PdfPages(pdf_filename) as pdf:
    fig = plt.figure(figsize=(5,10))
    for (patient_index, patient) in enumerate(patients):
        # Load patient data
        [samples, features, patient_timeseries] = read_patient_data(patient)
        nsamples = len(samples)
        nfeatures = len(features)
        projection = projected_data[patient_index]
        cmap = cm.rainbow(np.linspace(0,1,nsamples))

        ax = plt.subplot(subplot_ny, subplot_nx, patient_index+1)
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
            ax.plot(projection[0,0], projection[0,1], 'v', markersize=5)
            for i in range(nsamples-1):
                ax.plot(projection[i:i+2,0], projection[i:i+2,1], 'o-', markersize=3, linewidth=0.5)

            # Title subplot
            plt.title('patient %s' % patient, fontsize=9)

            # Fix up axes.
            plt.ylabel('tIC 2', fontsize=7)
            plt.xlabel('tIC 1', fontsize=7)
            plt.xticks([])
            plt.yticks([])
            plt.axis([-3, +3, -3, +3])

            plt.axis('square')

        except Exception as e:
            print(e)

    # Write figure
    pdf.savefig()
    plt.close()
