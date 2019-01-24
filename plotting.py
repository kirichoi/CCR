# -*- coding: utf-8 -*-
"""
CCR plotting module

Kiri Choi (c) 2018
"""

import os, sys
import tellurium as te
import roadrunner
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb


def plotAllProgress(listOfDistances, labels=None, SAVE_PATH=None):
    """
    Plots multiple convergence progress 
    
    :param listOfDistances: 2D array of distances
    :param labels: list of strings to use as labels
    :param SAVE: flag for saving the output
    """
    
    for i in range(len(listOfDistances)):
        plt.plot(listOfDistances[i])
    if labels:
        plt.legend(labels)
    plt.xlabel("Generations", fontsize=15)
    plt.ylabel("Distance", fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if SAVE_PATH is not None:
        if os.path.splitext(SAVE_PATH)[1] == '':
            plt.savefig(os.path.join(SAVE_PATH, 'images/convergence.pdf'), bbox_inches='tight')
        else:
            plt.savefig(SAVE_PATH, bbox_inches='tight')
    plt.show()

def plotProgress(distance, SAVE_PATH=None):
    """
    Plots convergence progress
    
    :param distance: array of distances
    :param model_type: reference model type, e.g. 'FFL', 'Linear', etc.
    :param SAVE: flag for saving the output
    """
    
    plt.plot(distance)
    plt.xlabel("Generations", fontsize=15)
    plt.ylabel("Distance", fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if SAVE_PATH is not None:
        if os.path.splitext(SAVE_PATH)[1] == '':
            plt.savefig(os.path.join(SAVE_PATH, 'images/convergence.pdf'), bbox_inches='tight')
        else:
            plt.savefig(SAVE_PATH, bbox_inches='tight')
    plt.show()

def plotResidual(realModel, ens_model, ens_dist, SAVE_PATH=None):
    """
    Plots residuals
    
    :param realModel: reference model
    :param ens_model: model ensemble
    :param ens_dist: model distances
    :param model_typ: reference model type
    :param SAVE: flag for saving the output
    """
    
    r_real = te.loada(realModel)
    result_real = r_real.simulate(0, 100, 100)
    
    top_result = []
    top_diff = []
    
    for i in range(len(ens_model)):
        r = te.loada(ens_model[np.argsort(ens_dist)[i]])
        top_sim = r.simulate(0, 100, 100)
        top_result.append(top_sim)
        top_diff.append(np.subtract(result_real[:,1:], top_sim[:,1:]))

    percentage = 0.1#float(pass_size)/ens_size
    
    ave_diff = np.average(top_diff[:int(len(ens_model)*percentage)], axis=0)
    
    plt.plot(ave_diff)
    plt.xlabel("Time (s)", fontsize=15)
    plt.ylabel("Residual", fontsize=15)
    plt.legend(r.getFloatingSpeciesIds())
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if SAVE_PATH is not None:
        if os.path.splitext(SAVE_PATH)[1] == '':
            plt.savefig(os.path.join(SAVE_PATH, 'images/average_residual.pdf'), bbox_inches='tight')
        else:
            plt.savefig(SAVE_PATH, bbox_inches='tight')
    plt.show()
    
def plotHistogram():
    """
    """
    
def plotDistanceHistogram(ens_dist, nbin=25, SAVE_PATH=None):
    """
    """
    
    plt.hist(ens_dist, bins=nbin, density=True)
    plt.xlabel("Distance", fontsize=15)
    plt.ylabel("Normalized Frequency", fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if SAVE_PATH is not None:
        if os.path.splitext(SAVE_PATH)[1] == '':
            plt.savefig(os.path.join(SAVE_PATH, 'images/distance_hist.pdf'), bbox_inches='tight')
        else:
            plt.savefig(SAVE_PATH, bbox_inches='tight')
    plt.show()

def plotDistanceHistogramWithKDE(dist_top, log_dens, minInd, nbin=40, SAVE_PATH=None):
    """
    """
    
    hist = plt.hist(dist_top, bins=nbin, density=True)
    plt.vlines(dist_top[minInd[0][0]], 0, np.max(hist[0]), linestyles='dashed')
    plt.plot(dist_top, np.exp(log_dens), color='tab:red')
    plt.xlabel("Distance", fontsize=15)
    plt.ylabel("Normalized Frequency", fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if SAVE_PATH is not None:
        if os.path.splitext(SAVE_PATH)[1] == '':
            plt.savefig(os.path.join(SAVE_PATH, 'images/distance_hist_w_KDE.pdf'), bbox_inches='tight')
        else:
            plt.savefig(SAVE_PATH, bbox_inches='tight')
    plt.show()

