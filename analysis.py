# -*- coding: utf-8 -*-
"""
CCR analysis extension module

Kiri Choi (c) 2018
"""

import tellurium as te
import roadrunner
import numpy as np
from scipy import signal
from sklearn import neighbors
import plotting as pt

def selectWithCutoff(model_top, dist_top, cutoff=0.1):
    """
    Model selection routine that returns a list of models with distances within
    the defined percentage.
    
    :param model_top: list of models sorted according to corresponding distances
    :param dist_top: list of sorted distances
    :param cutoff: percentage to cutoff
    """
    
    coind = int(len(model_top)*cutoff)
    pt.plot_distance_histogram(dist_top, cutoff_val=dist_top[coind])
    
    return model_top[:coind], dist_top[:coind]

def selectWithKernalDensity(model_top, dist_top):
    """
    Model selection rountine that returns a list of models based on the output
    of kernal density estimation.
    
    :param model_top: list of models sorted according to corresponding distances
    :param dist_top: list of sorted distances
    """
    
    dist_top_reshape = dist_top.reshape((100,1))
    
    kde = neighbors.KernelDensity(kernel='tophat', bandwidth=0.1).fit(dist_top_reshape)
    
    log_dens = kde.score_samples(dist_top_reshape)
    
    minInd = signal.argrelextrema(log_dens, np.less)
    
    return minInd, log_dens
