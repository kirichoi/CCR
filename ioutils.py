# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:18:02 2018

@author: ckiri
"""

import sys, os
import numpy as np
import pandas as pd

def exportSettings(settingsDict):
    """
    """
    
def exportOutputs(models, dists, best_dist, avg_dist, settings, time, ens_st, path=None):
    """
    """
    
    if path:
        outputdir = path
    else:
        outputdir = os.path.join(os.getcwd(), 'output')
        
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    if not os.path.exists(os.path.join(outputdir, 'models')):
        os.mkdir(os.path.join(outputdir, 'models'))
    
    df = pd.DataFrame(np.array(dists))
    df.to_csv(os.path.join(outputdir, 'dist_collected.txt'))
    
    stat = pd.DataFrame(np.array([best_dist, avg_dist]).T, columns=['generation best', 'generation average'])
    stat.to_csv(os.path.join(outputdir, 'dist_stat.txt'))
    
    outputtxt = open(os.path.join(outputdir, 'report.txt'), 'w')
    outputtxt.writelines('------------------------- REPORT -------------------------\n')
    outputtxt.writelines('RUN COMPLETE. HERE ARE SOME METRIC YOU MIGHT BE INTERESTED\n')
    outputtxt.writelines('No. of Generations: ' + str(settings['n_gen']) + '\n')
    outputtxt.writelines('Ensemble Size: ' + str(settings['ens_size']) + '\n')
    outputtxt.writelines('No. of Collected Models: ' + str(len(models)) + '\n')
    outputtxt.writelines('Run Time: ' + str(time) + ' s\n')
    outputtxt.writelines('No. Stoich. Analyzed: ' + str(len(ens_st)) + '\n')
    outputtxt.close()
    
    for i in range(len(models)):
        modeltxt = open(os.path.join(outputdir, 'models/model_' + str(i) + '.txt'), 'w')
        modeltxt.write(models[i])
        modeltxt.close()
    

def readSettings(settingsPath):
    """
    """


def readModels(modelsPath):
    """
    """
    
    modelfiles = [f for f in os.listdir(modelsPath) if os.path.isfile(os.path.join(modelsPath, f))]

    antstr = []
    for i in modelfiles:
        sbmlstr = open(os.path.join(modelsPath, i), 'r')
        antstr.append(sbmlstr.read())
        sbmlstr.close()
    
    return antstr


def readData(dataPath):
    """
    """
    
    
def testModels(modelType):
    """
    """
    
    if modelType == 'Linear':
        # Linear    
        realModel = """
        var S1, S2, S3, S4;
        const S0, S5;
        J0: S0 -> S1; k0*S0;
        J1: S1 -> S2; k1*S1;
        J2: S2 -> S3; k2*S2;
        J3: S3 -> S4; k3*S3;
        J4: S4 -> S5; k4*S4;
       
        k0 = 0.285822003905
        k1 = 0.571954691013
        k2 = 0.393173236422
        k3 = 0.75830845241
        k4 = 0.27503984992
       
        S0 = 4
        S5 = 5
        S1 = 0
        S2 = 0
        S3 = 0
        S4 = 0
        """
    elif modelType == 'Nested':
        # Nested
        realModel = """
        var S1, S2, S3;
        const S0, S4;
        J0: S0 -> S1; k0*S0;
        J1: S1 -> S2; k1*S1;
        J2: S2 -> S1; k2*S2;
        J3: S1 -> S3; k3*S1;
        J4: S3 -> S1; k4*S3;
        J5: S3 -> S4; k5*S3;
       
        k0 = 0.285822003905
        k1 = 0.571954691013
        k2 = 0.393173236422
        k3 = 0.75830845241
        k4 = 0.148522702962
        k5 = 0.348927696783
       
        S0 = 3
        S4 = 5
        S1 = 0
        S2 = 0
        S3 = 0
        """
    elif modelType == 'FFL':
        # FFL    
        realModel = """
        var S1, S2, S3;
        const S0, S4;
        J0: S0 -> S1; k0*S0;
        J1: S1 -> S2; k1*S1;
        J2: S2 -> S3; k2*S2;
        J3: S3 -> S4; k3*S3;
        J4: S1 -> S3; k4*S1;
       
        k0 = 0.285822003905
        k1 = 0.571954691013
        k2 = 0.393173236422
        k3 = 0.75830845241
        k4 = 0.148522702962
       
        S0 = 3
        S4 = 5
        S1 = 0
        S2 = 0
        S3 = 0
        """
    elif modelType == 'Branched':
        #Branched
        realModel = """
        var S1, S2, S3, S4, S5;
        const S0, S6;
        J0: S0 -> S1; k0*S0;
        J1: S1 -> S2; k1*S1;
        J2: S1 -> S3; k2*S1;
        J3: S3 -> S4; k3*S3;
        J4: S3 -> S5; k4*S3;
        J5: S2 -> S6; k5*S2;
        J6: S4 -> S6; k6*S4;
        J7: S5 -> S6; k7*S5;
       
        k0 = 0.285822003905
        k1 = 0.571954691013
        k2 = 0.393173236422
        k3 = 0.75830845241
        k4 = 0.148522702962
        k5 = 0.348927696783
        k6 = 0.572677236248
        k7 = 0.497208763889
       
        S0 = 4
        S6 = 3
        S1 = 0
        S2 = 0
        S3 = 0
        S4 = 0
        S5 = 0
        """
    else:
        raise Exception("Requested model not found")
        
    return realModel