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

def plot_progress(distance, model_type, SAVE=False):
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
    if SAVE:
        plt.savefig(os.path.join('./convergence_' + model_type + '.pdf'), bbox_inches='tight')
    plt.show()

def plot_residual(realModel, ens_model, ens_dist, model_type, SAVE=False):
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
    if SAVE:
        plt.savefig(os.path.join('./average_residual_' + model_type + '.pdf'), bbox_inches='tight')
    plt.show()
    
def plot_histogram():
    """
    """
    
def plot_distance_histogram(ens_dist, nbin=25, cutoff_val=None):
    """
    """
    
    plt.hist(ens_dist, bins=nbin, density=True)
    if cutoff_val is not None:
        plt.vlines(cutoff_val, )
    plt.xlabel("Distance", fontsize=15)
    plt.ylabel("Normalized Frequency", fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()
    

