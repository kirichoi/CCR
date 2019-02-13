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
    
def exportOutputs(models, dists, best_dist, avg_dist, settings, time, rl_track, path=None):
    """
    Export all outputs to specified path
    
    """
    
    if path:
        outputdir = path
    else:
        outputdir = os.path.join(os.getcwd(), 'output')
        
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    if not os.path.exists(os.path.join(outputdir, 'models')):
        os.mkdir(os.path.join(outputdir, 'models'))
    
    df = pd.DataFrame(np.array(dists), columns=['distance'])
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
    outputtxt.writelines('No. Stoich. Analyzed: ' + str(len(rl_track)) + '\n')
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
    Read model files
    
    :param modelsPath: path to a directory containing model files
    :returns: list of model strings
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
    Read data encoded in csv.
    
    :param dataPath: path to a csv file
    :returns: DataFrame
    """
    
    if os.path.exists(dataPath):
        df = pd.read_csv(dataPath)
        return df
    else:
        raise Exception("Cannot find the file at the specified path")
    
    
def testModels(modelType):
    """
    Returns a test model
    
    :param modelType: model name, e.g. 'FFL', 'Linear', 'Nested', 'Branched'
    :returns: Antimony string
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
    elif modelType == 'FFL_m':
        # FFL    
        realModel = """
        var S1, S2, S3;
        const S0, S4;
        J0: S0 -> S1; Kf0*S0/(1 + S0);
        J1: S1 -> S2; Kf1*S1/(1 + S1);
        J2: S2 -> S3; Kf2*S2/(1 + S2);
        J3: S3 -> S4; Kf3*S3/(1 + S3);
        J4: S1 -> S3; Kf4*S1/(1 + S1);
       
        Kf0 = 0.285822003905
        Kf1 = 0.571954691013
        Kf2 = 0.393173236422
        Kf3 = 0.75830845241
        Kf4 = 0.148522702962
       
        S0 = 3
        S4 = 5
        S1 = 0
        S2 = 0
        S3 = 0
        """
    elif modelType == 'FFL_r':
        # FFL    
        realModel = """
        var S1, S2, S3;
        const S0, S4;
        J0: S0 -> S1; (Kf0*S0 - Kr0*S1)/(1 + S0 + S1);
        J1: S1 -> S2; (Kf0*S1 - Kr1*S2)/(1 + S1 + S2);
        J2: S2 -> S3; (Kf2*S2 - Kr2*S3)/(1 + S2 + S3);
        J3: S3 -> S4; (Kf3*S3 - Kr3*S4)/(1 + S3 + S4);
        J4: S1 -> S3; (Kf4*S1 - Kr4*S3)/(1 + S1 + S3);
       
        Kf0 = 0.285822003905
        Kf1 = 0.571954691013
        Kf2 = 0.393173236422
        Kf3 = 0.75830845241
        Kf4 = 0.548522702962
        Kr0 = 0.072790656829
        Kr1 = 0.27828563882
        Kr2 = 0.166906190765
        Kr3 = 0.27057924345
        Kr4 = 0.1026847763
       
        S0 = 3
        S4 = 5
        S1 = 0
        S2 = 0
        S3 = 0
        """
    elif modelType == 'FFL_a':
        # FFL    
        realModel = """
        var S1, S2, S3;
        const S0, S4;
        J0: S0 -> S1; Kf0*S0/(1 + S0);
        J1: S1 -> S2; Kf1*S1/(1 + S1);
        J2: S2 -> S3; Kf2*S2/(1 + S2);
        J3: S3 -> S4; Kf3*S3/(1 + S3);
        J4: S1 -> S3; (Kf4*S1/(1 + S1))*(1 + Ka0*S2);
       
        Kf0 = 0.285822003905
        Kf1 = 0.571954691013
        Kf2 = 0.393173236422
        Kf3 = 0.75830845241
        Kf4 = 0.148522702962
        Ka0 = 0.883848629231
       
        S0 = 3
        S4 = 5
        S1 = 0
        S2 = 0
        S3 = 0
        """
    elif modelType == 'FFL_i':
        # I1FFL
        realModel = """
        var S1, S2, S3;
        const S0, S4;
        J0: S0 -> S1; Kf0*S0/(1 + S0);
        J1: S1 -> S2; Kf1*S1/(1 + S1);
        J2: S2 -> S3; Kf2*S2/(1 + S2);
        J3: S3 -> S4; Kf3*S3/(1 + S3);
        J4: S1 -> S3; (Kf4*S1/(1 + S1))*1/(1 + Ki0*S2);
       
        Kf0 = 0.285822003905
        Kf1 = 0.571954691013
        Kf2 = 0.393173236422
        Kf3 = 0.75830845241
        Kf4 = 0.148522702962
        Ki0 = 0.974569278466
       
        S0 = 3
        S4 = 5
        S1 = 0
        S2 = 0
        S3 = 0
        """
    elif modelType == 'FFL_r_a':
        # FFL    
        realModel = """
        var S1, S2, S3;
        const S0, S4;
        J0: S0 -> S1; (Kf0*S0 - Kr0*S1)/(1 + S0 + S1);
        J1: S1 -> S2; (Kf0*S1 - Kr1*S2)/(1 + S1 + S2);
        J2: S2 -> S3; (Kf2*S2 - Kr2*S3)/(1 + S2 + S3);
        J3: S3 -> S4; (Kf3*S3 - Kr3*S4)/(1 + S3 + S4);
        J4: S1 -> S3; (1 + Ka0*S2)*(Kf4*S1 - Kr4*S3)/(1 + S1 + S3);
       
        Kf0 = 0.285822003905
        Kf1 = 0.571954691013
        Kf2 = 0.393173236422
        Kf3 = 0.75830845241
        Kf4 = 0.548522702962
        Kr0 = 0.072790656829
        Kr1 = 0.27828563882
        Kr2 = 0.166906190765
        Kr3 = 0.27057924345
        Kr4 = 0.1026847763
        Ka0 = 0.883848629231
       
        S0 = 3
        S4 = 5
        S1 = 0
        S2 = 0
        S3 = 0
        """
    elif modelType == 'FFL_r_i':
        # FFL    
        realModel = """
        var S1, S2, S3;
        const S0, S4;
        J0: S0 -> S1; (Kf0*S0 - Kr0*S1)/(1 + S0 + S1);
        J1: S1 -> S2; (Kf0*S1 - Kr1*S2)/(1 + S1 + S2);
        J2: S2 -> S3; (Kf2*S2 - Kr2*S3)/(1 + S2 + S3);
        J3: S3 -> S4; (Kf3*S3 - Kr3*S4)/(1 + S3 + S4);
        J4: S1 -> S3; 1/(1 + Ki0*S2)*(Kf4*S1 - Kr4*S3)/(1 + S1 + S3);
       
        Kf0 = 0.285822003905
        Kf1 = 0.571954691013
        Kf2 = 0.393173236422
        Kf3 = 0.75830845241
        Kf4 = 0.548522702962
        Kr0 = 0.072790656829
        Kr1 = 0.27828563882
        Kr2 = 0.166906190765
        Kr3 = 0.27057924345
        Kr4 = 0.1026847763
        Ki0 = 0.974569278466
       
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
