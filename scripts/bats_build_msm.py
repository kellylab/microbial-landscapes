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

filename = '../data/bats_orig.txt'

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

nstates = 5 # number of hidden states

# for depth in depths:
#     print(depth)
#     features = np.log(np.array(df[df['depth'] == depth][feature_names])) # log populations
#
#     X_sklearn = features  # sklearn style input: (n_samples, n_features)
#     X_msmb = [X_sklearn]  # MSMBuilder style input: list of (n_samples, n_features)
#
#
#     clusterer_sklearn = sklearn.cluster.KMeans(n_clusters=nstates)
#     clusterer_sklearn.fit(X_sklearn)
#
#     clusterer_msmb = msmbuilder.cluster.KMeans(n_clusters=nstates)
#     clusterer_msmb.fit(X_msmb)
#
#     ghmm = msmbuilder.hmm.GaussianHMM(n_states=nstates)
#     ghmm.fit(X_msmb)

#
# MAKE PLOTS
#

# Get colormap
from matplotlib.pyplot import cm
cmap = cm.rainbow(np.linspace(0,1,ndepths))

# Make plots
pdf_filename = '../figures/bats_tica.pdf'
with PdfPages(pdf_filename) as pdf:
    # Plot individual trajectories
    plt.figure(figsize=(10,10))

    for (subplot_index, feature_name) in enumerate(feature_names):
        plt.subplot(nfeatures, 1, subplot_index+1)

        plt.hold(True)
        for (depth, color) in zip(depths, cmap):
            features = np.log(np.array(df[df['depth'] == depth][feature_names])) # log abundances
            plt.plot(features[:,subplot_index], '-', color=color)

        plt.title(feature_names[subplot_index])
        plt.xlabel('measurement number')
        if subplot_index+1 == nfeatures: plt.ylabel('log abundance')
    plt.legend(depths)
    pdf.savefig()
    plt.close()

    # Plot joint trajectories
    plt.figure(figsize=(10,10))
    plt.hold(True)
    ax = plt.gca()

    ellipses = list()
    for (depth, color) in zip(depths, cmap):
        print(depth)

        features = np.log(np.array(df[df['depth'] == depth][feature_names]))
        plt.plot(features[:,0], features[:,1], '-', color=color, alpha=0.3)

        try:
            ghmm = msmbuilder.hmm.GaussianHMM(n_states=nstates, n_iter=500, n_lqa_iter=50, reversible_type='mle')
            ghmm.fit([features])
            #plt.plot(ghmm.means_[:,0], ghmm.means_[:,1], 'o', color=color)
            print ghmm.means_
            print ghmm.vars_
            from matplotlib.patches import Ellipse
            for state in range(nstates):
                ellipse = Ellipse((ghmm.means_[state,0], ghmm.means_[state,1]), width=np.sqrt(ghmm.vars_[state,0]), height=np.sqrt(ghmm.vars_[state,1]), color=color)
                ellipse.set_alpha(1.0)
                ellipses.append(ellipse)
                #ax.add_patch(ellipse)
        except Exception as e:
            print(e)

    ax = plt.gca()
    for ellipse in ellipses:
        ax.add_patch(ellipse)

    plt.xlabel(feature_names[0])
    plt.ylabel(feature_names[1])
    plt.legend(depths)
    pdf.savefig()
    plt.close()

    # Plot 3D
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.hold(True)
    for (depth, color) in zip(depths, cmap):
        features = np.log(np.array(df[df['depth'] == depth][feature_names]))
        ax.plot(features[:,0], features[:,1], features[:,2], '-', color=color)

    ax.set_xlabel(feature_names[0])
    ax.set_ylabel(feature_names[1])
    ax.set_zlabel(feature_names[2])
    ax.legend(depths)
    pdf.savefig()
    plt.close()
