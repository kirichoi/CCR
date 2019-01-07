# -*- coding: utf-8 -*-
"""
Created on Sat May 19 16:01:19 2018

@author: kirichoi
"""

import os, sys
import tellurium as te
import roadrunner
import numpy as np
import antimony
import scipy.optimize
import networkGenerator as ng
import plotting as pt
import ioutils
import analysis
import matplotlib.pyplot as plt
import time

#np.seterr(all='raise')

def f1(k_list, *args):
    global counts
    global countf
    
    args[0].reset()
    
    args[0].setValues(args[0].getGlobalParameterIds(), k_list)
    
#    ss = args[0].steadyStateSolver
#    ss.allow_approx = True
#    ss.allow_presimulation = False
    
    try:
        r.steadyState()
        objCC = args[0].getScaledConcentrationControlCoefficientMatrix()
        dist_obj = (np.linalg.norm(realConcCC - objCC))
    except:
        countf += 1
        dist_obj = 10000
        
    counts += 1
    
    return dist_obj


def callbackF(X, convergence=0.):
    global counts
    global countf
    print(str(counts) + ", " + str(countf))
    return False


def mutate_and_evaluate(listantStr, listdist):
    global countf
    global counts
    
    eval_dist = np.empty(mut_size)
    eval_model = np.empty(mut_size, dtype='object')
    
    for m in mut_range:
        antimony.loadAntimonyString(listantStr[m])
        module = antimony.getModuleNames()[-1]
        
        r = te.loada(listantStr[m])
        param_val = r.getGlobalParameterValues()
        r.steadyState()
        
        tempdiff = np.max(np.abs(realConcCC - 
                r.getScaledConcentrationControlCoefficientMatrix()), axis=0)
        
        stt = [[],[],[]]
        st = ens_st[0]
        
        o = 0
        
        while ((stt[1] != realFloatingIdsInd or stt[2] != realBoundaryIdsInd or
                st in ens_st) and (o < maxIter_mut)):
            rct = np.array(antimony.getReactantNames(module)).tolist()
            prd = np.array(antimony.getProductNames(module)).tolist()
            
            r_idx = np.random.choice(np.arange(nr), p=np.divide(tempdiff,np.sum(tempdiff)))
            rt = ng.pickReactionType()
            if rt == ng.TReactionType.UNIUNI:
                # TODO: pick reactants and products based on control coefficients
                # UniUni
                posRctInd = np.append(np.array(realFloatingIdsInd)[np.nonzero(
                        realConcCC[:,r_idx])[0]], np.array(realBoundaryIdsInd))
                rct_id = np.random.choice(posRctInd, size=1)
                prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                
                # Search for potentially identical reactions
                all_rct = [i for i, x in enumerate(rct) if x == ['S'+str(rct_id[0])]]
                all_prd = [i for i, x in enumerate(prd) if x == ['S'+str(prd_id[0])]]
                
                while (len(set(all_rct) & set(all_prd)) > 0):
                    rct_id = np.random.choice(np.arange(ns), size=1)
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                    all_rct = [i for i, x in enumerate(rct) if x == ['S'+str(rct_id[0])]]
                    all_prd = [i for i, x in enumerate(prd) if x == ['S'+str(prd_id[0])]]
                
                rct[r_idx] = ["S" + str(rct_id[0])]
                prd[r_idx] = ["S" + str(prd_id[0])]
                
            elif rt == ng.TReactionType.BIUNI:
                # BiUni
                posRctInd = np.append(np.array(realFloatingIdsInd)[np.nonzero(
                        realConcCC[:,r_idx])[0]], np.array(realBoundaryIdsInd))
                # Pick two reactants
                rct_id = np.random.choice(posRctInd, size=2, replace=True)
                # pick a product but only products that don't include the reactants
                prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                
                 # Search for potentially identical reactions
                all_rct = [i for i, x in enumerate(rct) if x == ('S'+str(rct_id[0]),
                                    'S'+str(rct_id[1])) or x == ('S'+str(rct_id[1]), 
                                    'S'+str(rct_id[0]))]
                all_prd = [i for i, x in enumerate(prd) if x == ['S'+str(prd_id[0])]]
                
                while (len(set(all_rct) & set(all_prd)) > 0):
                    rct_id = np.random.choice(np.arange(ns), size=2, replace=True)
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                    all_rct = [i for i, x in enumerate(rct) if x == ('S'+str(rct_id[0]),
                                        'S'+str(rct_id[1])) or x == ('S'+str(rct_id[1]), 
                                        'S'+str(rct_id[0]))]
                    all_prd = [i for i, x in enumerate(prd) if x == ['S'+str(prd_id[0])]]
                
                rct[r_idx] = ["S" + str(rct_id[0]), "S" + str(rct_id[1])]
                prd[r_idx] = ["S" + str(prd_id[0])]
                
            elif rt == ng.TReactionType.UNIBI:
                # UniBi
                posRctInd = np.append(np.array(realFloatingIdsInd)[np.nonzero(
                        realConcCC[:,r_idx])[0]], np.array(realBoundaryIdsInd))
                rct_id = np.random.choice(posRctInd, size=1)
                # pick a product but only products that don't include the reactant
                prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2,
                                          replace=True)
                
                # Search for potentially identical reactions
                all_rct = [i for i, x in enumerate(rct) if x == ['S'+str(rct_id[0])]]
                all_prd = [i for i, x in enumerate(prd) if x == ('S'+str(prd_id[0]),
                                    'S'+str(prd_id[1])) or x == ('S'+str(prd_id[1]), 
                                    'S'+str(prd_id[0]))]
                
                while (len(set(all_rct) & set(all_prd)) > 0):
                    rct_id = np.random.choice(np.arange(ns), size=1)
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), 
                                              size=2, replace=True)
                    all_rct = [i for i, x in enumerate(rct) if x == ['S'+str(rct_id[0])]]
                    all_prd = [i for i, x in enumerate(prd) if x == ('S'+str(prd_id[0]),
                                        'S'+str(prd_id[1])) or x == ('S'+str(prd_id[1]), 
                                        'S'+str(prd_id[0]))]
                
                rct[r_idx] = ["S" + str(rct_id[0])]
                prd[r_idx] = ["S" + str(prd_id[0]), "S" + str(prd_id[1])]
                
            else:
                # BiBi
                posRctInd = np.append(np.array(realFloatingIdsInd)[np.nonzero(
                        realConcCC[:,r_idx])[0]], np.array(realBoundaryIdsInd))
                rct_id = np.random.choice(posRctInd, size=2, replace=True)
                # pick a product but only products that don't include the reactant
                prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, 
                                          replace=True)
                
                # Search for potentially identical reactions
                all_rct = [i for i, x in enumerate(rct) if x == ('S'+str(rct_id[0]),
                                    'S'+str(rct_id[1])) or x == ('S'+str(rct_id[1]),
                                    'S'+str(rct_id[0]))]
                all_prd = [i for i, x in enumerate(prd) if x == ('S'+str(prd_id[0]),
                                    'S'+str(prd_id[1])) or x == ('S'+str(prd_id[1]), 
                                    'S'+str(prd_id[0]))]
                
                while (len(set(all_rct) & set(all_prd)) > 0):
                    rct_id = np.random.choice(np.arange(ns), size=2, replace=True)
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id),
                                              size=2, replace=True)
                    all_rct = [i for i, x in enumerate(rct) if x == ('S'+str(rct_id[0]),
                                        'S'+str(rct_id[1])) or x == ('S'+str(rct_id[1]),
                                        'S'+str(rct_id[0]))]
                    all_prd = [i for i, x in enumerate(prd) if x == ('S'+str(prd_id[0]),
                                        'S'+str(prd_id[1])) or x == ('S'+str(prd_id[1]),
                                        'S'+str(prd_id[0]))]
                
                rct[r_idx] = ["S" + str(rct_id[0]), "S" + str(rct_id[1])]
                prd[r_idx] = ["S" + str(prd_id[0]), "S" + str(prd_id[1])]
                
            reactionList = []
            for i in r_range:
                if len(rct[i]) == 1 and len(prd[i]) == 1:
                    rtype = ng.TReactionType.UNIUNI
                elif len(rct[i]) == 1 and len(prd[i]) == 2:
                    rtype = ng.TReactionType.UNIBI
                elif len(rct[i]) == 2 and len(prd[i]) == 1:
                    rtype = ng.TReactionType.BIUNI
                else:
                    rtype = ng.TReactionType.BIBI
                
                rct_ind = [s.replace('S', '') for s in rct[i]]
                rct_ind = list(map(int, rct_ind))
                prd_ind = [s.replace('S', '') for s in prd[i]]
                prd_ind = list(map(int, prd_ind))
                reactionList.append([rtype, rct_ind, prd_ind, param_val[i]]) 
            
            st = ng.getFullStoichiometryMatrix(reactionList, ns).tolist()
            stt = ng.removeBoundaryNodes(np.array(st))
            o += 1
        
        if o >= maxIter_mut:
            eval_dist[m] = listdist[m]
            eval_model[m] = listantStr[m]
        else:
            antStr = ng.generateAntimony(realFloatingIds, realBoundaryIds, 
                                          stt[1], stt[2], reactionList, 
                                          boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)
                
                counts = 0
                countf = 0
                
                ss = r.steadyStateSolver
                ss.allow_approx = True
                ss.allow_presimulation = False
                r.steadyState()
                
                p_bound = [(1e-3, 1.)]*r.getNumGlobalParameters()
                res = scipy.optimize.differential_evolution(f1, args=(r,), 
                            bounds=p_bound, maxiter=optiMaxIter, tol=optiTol,
                            polish=optiPolish, seed=r_seed)
                if not res.success:
                    eval_dist[m] = listdist[m]
                    eval_model[m] = listantStr[m]
                else:
                    r = te.loada(antStr)
                    r.setValues(r.getGlobalParameterIds(), res.x)
                    
                    r.steadyState()
                    SS_i = r.getFloatingSpeciesConcentrations()
                    # Buggy model
                    if np.any(SS_i > 1e5):
                        r.reset()
                        ss.allow_presimulation = True
                        ss.presimulation_time = 100
                        r.steadyState()
                        SS_i = r.getFloatingSpeciesConcentrations()
                    if np.any(SS_i < 1e-5) or np.any(SS_i > 1e5):
                        eval_dist[m] = listdist[m]
                        eval_model[m] = listantStr[m]
                    else:
                        concCC_i = r.getScaledConcentrationControlCoefficientMatrix()
                        
                        if np.isnan(concCC_i).any():
                            eval_dist[m] = listdist[m]
                            eval_model[m] = listantStr[m]
                        else:
                            concCC_i_row = concCC_i.rownames
                            concCC_i_col = concCC_i.colnames
                            concCC_i = concCC_i[np.argsort(concCC_i_row)]
                            concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                            
                            dist_i = w1*(np.linalg.norm(realConcCC - concCC_i))
                            
                            if dist_i < listdist[m]:
                                eval_dist[m] = dist_i
                                r.reset()
                                eval_model[m] = r.getAntimony(current=True)
                                ens_st.append(st)
                            else:
                                eval_dist[m] = listdist[m]
                                eval_model[m] = listantStr[m]
            except:
                eval_dist[m] = listdist[m]
                eval_model[m] = listantStr[m]
        antimony.clearPreviousLoads()

    return eval_dist, eval_model


def initialize():
    global countf
    global counts
    
    numBadModels = 0
    numGoodModels = 0
    numIter = 0
    
    ens_dist = np.empty(ens_size)
    ens_model = np.empty(ens_size, dtype='object')
    ens_st = []
    
    # Initial Random generation
    while (numGoodModels < ens_size):
        rl = ng.generateReactionList(ns, nr)
        st = ng.getFullStoichiometryMatrix(rl, ns).tolist()
        stt = ng.removeBoundaryNodes(np.array(st))
        # Ensure no redundant model
        while (stt[1] != realFloatingIdsInd or stt[2] != realBoundaryIdsInd or st in ens_st):
            rl = ng.generateReactionList(ns, nr)
            st = ng.getFullStoichiometryMatrix(rl, ns).tolist()
            stt = ng.removeBoundaryNodes(np.array(st))
        antStr = ng.generateAntimony(realFloatingIds, realBoundaryIds, stt[1],
                                      stt[2], rl, boundary_init=realBoundaryVal)
        try:
            r = te.loada(antStr)

            counts = 0
            countf = 0
            
            ss = r.steadyStateSolver
            ss.allow_approx = True
            ss.allow_presimulation = False
            r.steadyState()
            
            p_bound = [(1e-3, 1.)]*r.getNumGlobalParameters()
            res = scipy.optimize.differential_evolution(f1, args=(r,), 
                               bounds=p_bound, maxiter=optiMaxIter, tol=optiTol,
                               polish=optiPolish, seed=r_seed)
            if not res.success:
                numBadModels += 1
            else:
                # TODO: Might be able to cut the bottom part by simply using the obj func value from optimizer
                r = te.loada(antStr)
                r.setValues(r.getGlobalParameterIds(), res.x)
                    
                r.steadyState()
                SS_i = r.getFloatingSpeciesConcentrations()
                # Buggy model
                if np.any(SS_i > 1e5):
                    r.reset()
                    ss.allow_presimulation = True
                    ss.presimulation_time = 100
                    r.steadyState()
                    SS_i = r.getFloatingSpeciesConcentrations()
                        
                if np.any(SS_i < 1e-5) or np.any(SS_i > 1e5):
                    numBadModels += 1
                else:
                    concCC_i = r.getScaledConcentrationControlCoefficientMatrix()
                    
                    if np.isnan(concCC_i).any():
                        numBadModels += 1
                    else:
                        concCC_i_row = concCC_i.rownames
                        concCC_i_col = concCC_i.colnames
                        concCC_i = concCC_i[np.argsort(concCC_i_row)]
                        concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                        
                        dist_i = w1*(np.linalg.norm(realConcCC - concCC_i))
                        
                        ens_dist[numGoodModels] = dist_i
                        r.reset()
                        ens_model[numGoodModels] = r.getAntimony(current=True)
                        ens_st.append(st)
                        
                        numGoodModels = numGoodModels + 1
        except:
            numBadModels = numBadModels + 1
        antimony.clearPreviousLoads()
        numIter = numIter + 1
        if int(numIter/1000) == (numIter/1000):
            print("Number of iterations = " + str(numIter))
        if int(numIter/10000) == (numIter/10000):
            print("Number of good models = " + str(numGoodModels))
    
    print("In generation: 1")
    print("Number of total iterations = " + str(numIter))
    print("Number of bad models = " + str(numBadModels))
    
    return ens_dist, ens_model, ens_st


def random_gen(listAntStr, listDist):
    global countf
    global counts
    
    rndSize = len(listDist)
    
    rnd_dist = np.empty(rndSize)
    rnd_model = np.empty(rndSize, dtype='object')
    
    for l in range(rndSize):
        d = 0
        rl = ng.generateReactionList(ns, nr)
        st = ng.getFullStoichiometryMatrix(rl, ns).tolist()
        stt = ng.removeBoundaryNodes(np.array(st))
        # Ensure no redundant models
        while ((stt[1] != realFloatingIdsInd or stt[2] != realBoundaryIdsInd or
                st in ens_st) and (d < maxIter_gen)):
            rl = ng.generateReactionList(ns, nr)
            st = ng.getFullStoichiometryMatrix(rl, ns).tolist()
            stt = ng.removeBoundaryNodes(np.array(st))
            d += 1
        if d >= maxIter_gen:
            rnd_dist[l] = listDist[l]
            rnd_model[l] = listAntStr[l]
        else:
            antStr = ng.generateAntimony(realFloatingIds, realBoundaryIds, 
                            stt[1], stt[2], rl, boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)
                
                counts = 0
                countf = 0
                
                ss = r.steadyStateSolver
                ss.allow_approx = True
                ss.allow_presimulation = False
                r.steadyState()
                
                p_bound = [(1e-3, 1.)]*r.getNumGlobalParameters()
                res = scipy.optimize.differential_evolution(f1, args=(r,), 
                            bounds=p_bound, maxiter=optiMaxIter, tol=optiTol,
                            polish=optiPolish, seed=r_seed)
                # Failed to find solution
                if not res.success:
                    rnd_dist[l] = listDist[l]
                    rnd_model[l] = listAntStr[l]
                else:
                    r = te.loada(antStr)
                    r.setValues(r.getGlobalParameterIds(), res.x)
                    
                    r.steadyState()
                    SS_i = r.getFloatingSpeciesConcentrations()
                    # Buggy model
                    if np.any(SS_i > 1e5):
                        r.reset()
                        ss.allow_presimulation = True
                        ss.presimulation_time = 100
                        r.steadyState()
                        SS_i = r.getFloatingSpeciesConcentrations()
                    if np.any(SS_i < 1e-5) or np.any(SS_i > 1e5):
                        rnd_dist[l] = listDist[l]
                        rnd_model[l] = listAntStr[l]
                    else:
                        concCC_i = r.getScaledConcentrationControlCoefficientMatrix()
                        
                        if np.isnan(concCC_i).any():
                            rnd_dist[l] = listDist[l]
                            rnd_model[l] = listAntStr[l]
                        else:
                            concCC_i_row = concCC_i.rownames
                            concCC_i_col = concCC_i.colnames
                            concCC_i = concCC_i[np.argsort(concCC_i_row)]
                            concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                            
                            dist_i = w1*(np.linalg.norm(realConcCC - concCC_i))
                            
                            if dist_i < listDist[l]:
                                rnd_dist[l] = dist_i
                                r.reset()
                                rnd_model[l] = r.getAntimony(current=True)
                                ens_st.append(st)
                            else:
                                rnd_dist[l] = listDist[l]
                                rnd_model[l] = listAntStr[l]
            except:
                rnd_dist[l] = listDist[l]
                rnd_model[l] = listAntStr[l]
        antimony.clearPreviousLoads()
    
    return rnd_dist, rnd_model

# TODO: simulated annealing (multiply to fitness for rate constants)
if __name__ == '__main__':
    roadrunner.Config.setValue(roadrunner.Config.ROADRUNNER_DISABLE_WARNINGS, 3)


#%% Settings
    
    # Input data
    INPUT = None
    READ_SETTINGS = None
    
    # Test models
    modelType = 'FFL' # 'FFL', 'Linear', 'Nested', 'Branched'
    
    # General settings
    n_gen = 400 # Number of generations
    ens_size = 100 # Size of output ensemble
    pass_size = int(ens_size/10) # Number of models passed on the next generation without mutation
    mut_size = int(ens_size/2) # Number of models to mutate
    maxIter_gen = 5000 # Maximum iteration allowed for random generation
    maxIter_mut = 5000 # Maximum iteration allowed for mutation
    
    # Optimizer settings
    optiMaxIter = 1000 # Maximum iteraction allowed for optimizer
    optiTol = 0.01
    optiPolish = True
    w1 = 1.0 # Weight for control coefficients when calculating the distance
    w2 = 1.0 # Weight for steady-state and flux when calculating the distance
    
    # Random settings
    r_seed = 123123 # random seed
    r_roulette = False # Flag for using random roulette or best of pair for selection process
    NOISE = False # Flag for adding Gaussian noise to steady-state and control coefficiant values
    ABS_NOISE_STD = 0.01 # Standard deviation of Gaussian noise
    REL_NOISE_STD = 0.1 # Standard deviation of Gaussian noise
    
    # Plotting settings
    PLOT = True # Flag for plots
    SAVE_PLOT = True # Flag for saving plots
    
    # Data settings
    EXPORT_OUTPUT = True # Flag for saving collected models
    EXPORT_SETTINGS = False # Flag for saving current settings
    EXPORT_PATH = './output_ffl' # Path to save the output
    
    # Flag to run algorithm
    RUN = True
    
    #%%
    if RUN:
        # Using one of the test models
        realModel = ioutils.testModels(modelType)
        
        realRR = te.loada(realModel)
        
        realSS = realRR.steadyStateSolver
        realSS.allow_approx = False
        realSS.allow_presimulation = False
        realSS.presimulation_time = 100
        
        realNumBoundary = realRR.getNumBoundarySpecies()
        realNumFloating = realRR.getNumFloatingSpecies()
        realFloatingIds = np.sort(realRR.getFloatingSpeciesIds())
        realFloatingIdsInd = list(map(int, [s.strip('S') for s in realRR.getFloatingSpeciesIds()]))
        realBoundaryIds = np.sort(realRR.getBoundarySpeciesIds())
        realBoundaryIdsInd = list(map(int,[s.strip('S') for s in realRR.getBoundarySpeciesIds()]))
        realBoundaryVal = realRR.getBoundarySpeciesConcentrations()
        realGlobalParameterIds = realRR.getGlobalParameterIds()
        
        realRR.steadyState()
        realSteadyState = realRR.getFloatingSpeciesConcentrations()
        realSteadyStateRatio = np.divide(realSteadyState, np.min(realSteadyState))
        realFlux = realRR.getReactionRates()
        realRR.reset()
        realRR.steadyState()
        realFluxCC = realRR.getScaledFluxControlCoefficientMatrix()
        realConcCC = realRR.getScaledConcentrationControlCoefficientMatrix()
        
        # Ordering
        realFluxCCrow = realFluxCC.rownames
        realFluxCCcol = realFluxCC.colnames
        realFluxCC = realFluxCC[np.argsort(realFluxCCrow)]
        realFluxCC = realFluxCC[:,np.argsort(realFluxCCcol)]
        
        realConcCCrow = realConcCC.rownames
        realConcCCcol = realConcCC.colnames
        realConcCC = realConcCC[np.argsort(realConcCCrow)]
        realConcCC = realConcCC[:,np.argsort(realConcCCcol)]
        
        ns = realNumBoundary + realNumFloating # Number of species
        nr = realRR.getNumReactions() # Number of reactions
        
        print("Original Control Coefficients")
        print(realConcCC)
        print("Original Steady State Ratio")
        print(realSteadyStateRatio)
        
        if NOISE:
            for i in range(len(realConcCC)):
                for j in range(len(realConcCC[i])):
                    realConcCC[i][j] = (realConcCC[i][j] + np.random.normal(0,ABS_NOISE_STD) 
                    + np.random.normal(0,np.abs(realConcCC[i][j]*REL_NOISE_STD)))
            
            print("Control Coefficients with Noise Added")
            print(realConcCC)
    #%%
        # Define seed and ranges
        np.random.seed(r_seed)
        
        best_dist = []
        avg_dist = []
        
        n_range = range(1, n_gen)
        ens_range = range(ens_size)
        mut_range = range(mut_size)
        r_range = range(nr)
        
    #%%
        t1 = time.time()
        
        # Initialize
        ens_dist, ens_model, ens_st = initialize()
        
        dist_top_ind = np.argsort(ens_dist)
        dist_top = ens_dist[dist_top_ind]
        model_top = ens_model[dist_top_ind]
        
        print("Minimum distance: " + str(dist_top[0]))
        print("Average distance: " + str(np.average(dist_top)))
        best_dist.append(dist_top[0])
        avg_dist.append(np.average(dist_top))
        
        breakFlag = False
        
        # TODO: Remove for loop
        for n in n_range:
            if r_roulette:
                ens_inv = np.divide(1, dist_top[pass_size:])
                ens_prob = np.divide(ens_inv, np.sum(ens_inv))
                mut_ind = np.random.choice(np.arange(pass_size, ens_size), 
                                           size=mut_size, replace=False, p=ens_prob)
            else:
                minind = np.argsort(ens_dist)[:pass_size]
                tarind = np.delete(np.arange(ens_size), minind)
                mut_p = 1/ens_dist[tarind]/np.sum(1/ens_dist[tarind])
                mut_ind = np.random.choice(tarind, size=mut_size-pass_size, 
                                                   replace=False, p=mut_p)
                mut_ind = np.append(mut_ind, minind)
                mut_ind_inv = np.setdiff1d(np.arange(ens_size), mut_ind)
            
            eval_dist, eval_model = mutate_and_evaluate(ens_model[mut_ind], ens_dist[mut_ind])
            
            ens_model[mut_ind] = eval_model
            ens_dist[mut_ind] = eval_dist
            
    #        for tt in range(len(mut_ind)):
    #            r = te.loada(ens_model[mut_ind[tt]])
    #            ss = r.steadyStateSolver
    #            ss.allow_approx = True
    #            ss.allow_presimulation = False
    #            try:
    #                r.steadyState()
    #            except:
    #                print("Failure detacted at mutation: ", mut_ind[tt])
    #                print(np.sort(mut_ind))
    #                breakFlag = True
    #                break
    #        
    #        if breakFlag:
    #            break
            
            rnd_dist, rnd_model = random_gen(ens_model[mut_ind_inv], ens_dist[mut_ind_inv])
            ens_model[mut_ind_inv] = rnd_model
            ens_dist[mut_ind_inv] = rnd_dist
            
            dist_top_ind = np.argsort(ens_dist)
            dist_top = ens_dist[dist_top_ind]
            model_top = ens_model[dist_top_ind]
            
            print("In generation: " + str(n + 1))
            print("Minimum distance: " + str(dist_top[0]))
            print("Average distance: " + str(np.average(dist_top)))
            best_dist.append(dist_top[0])
            avg_dist.append(np.average(dist_top))
            
    #        for tt in range(len(mut_ind_inv)):
    #            r = te.loada(ens_model[mut_ind_inv[tt]])
    #            ss = r.steadyStateSolver
    #            ss.allow_approx = True
    #            ss.allow_presimulation = False
    #            try:
    #                r.steadyState()
    #            except:
    #                print("Failure detacted at random gen: ", mut_ind_inv[tt])
    #                print(np.sort(mut_ind_inv))
    #                breakFlag = True
    #                break
    #        
    #        if breakFlag:
    #            break
            
            # Error check
            #if np.average(dist_top) > 10000:
                #break
    
        # Check run time
        t2 = time.time()
        print(t2 - t1)
        
        #%%
        # Collect models
        minInd, log_dens = analysis.selectWithKernalDensity(model_top, dist_top)
        model_col = model_top[:minInd[0][0]]
        dist_col = dist_top[:minInd[0][0]]
            
    #%%
        EXPORT_PATH = os.path.abspath(os.path.join(os.getcwd(), EXPORT_PATH))
        
        if PLOT:
            # Convergence
            if SAVE_PLOT:
                if not os.path.exists(EXPORT_PATH):
                    os.mkdir(EXPORT_PATH)
                if not os.path.exists(os.path.join(EXPORT_PATH, 'images')):
                    os.mkdir(os.path.join(EXPORT_PATH, 'images'))
                pt.plotProgress(best_dist, SAVE_PATH=os.path.join(EXPORT_PATH, 'images/convergence_best.pdf'))
                pt.plotProgress(avg_dist, SAVE_PATH=os.path.join(EXPORT_PATH, 'images/convergence_avg.pdf'))
            else:
                pt.plotProgress(best_dist)
                pt.plotProgress(avg_dist)
            # TODO: Add polishing with fast optimizer 
            
            # Average residual
            if SAVE_PLOT:
                pt.plotResidual(realModel, ens_model, ens_dist, SAVE_PATH=os.path.join(EXPORT_PATH, 'images/average_residual.pdf'))
            else:
                pt.plotResidual(realModel, ens_model, ens_dist)
                
            # Distance histogram with KDE
            if SAVE_PLOT:
                pt.plotDistanceHistogramWithKDE(dist_top, log_dens, minInd, SAVE_PATH=os.path.join(EXPORT_PATH, 'images/distance_hist_w_KDE.pdf'))
            else:
                pt.plotDistanceHistogramWithKDE(dist_top, log_dens, minInd)
                
            # RMSE histogram
            r_real = te.loada(realModel)
            k_real = r_real.getGlobalParameterValues()
            
            top_result_k = []
            top_diff_k = []
            
            for i in ens_range:
                r = te.loada(ens_model[np.argsort(ens_dist)[i]])
                top_k = r.getGlobalParameterValues()
                top_result_k.append(top_k)
                try:
                    top_diff_k.append(np.sqrt(np.divide(np.sum(np.square(np.subtract(
                            k_real, top_k))),len(k_real))))
                except:
                    top_diff_k.append(np.sqrt(np.divide(np.sum(np.square(np.subtract(
                            k_real, top_k[1:]))),len(k_real))))
            
            krmse = top_diff_k[:pass_size]
            
            plt.hist(krmse, bins=15, density=True)
            plt.xlabel("RMSE", fontsize=15)
            plt.ylabel("Normalized Frequency", fontsize=15)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            if SAVE_PLOT:
                plt.savefig(os.path.join(EXPORT_PATH, 'images/parameter_rmse_.pdf'), bbox_inches='tight')
            plt.show()
            
    #%%
        if EXPORT_SETTINGS or EXPORT_OUTPUT:
            settings = {}
            settings['n_gen'] = n_gen
            settings['ens_size'] = ens_size
            settings['pass_size'] = pass_size
            settings['mut_size'] = mut_size
            settings['maxIter_gen'] = maxIter_gen
            settings['maxIter_mut'] = maxIter_mut
            settings['optiMaxIter'] = optiMaxIter
            settings['optiTol'] = optiTol
            settings['optiPolish'] = optiPolish
            settings['r_seed'] = r_seed
            
            if EXPORT_SETTINGS:
                ioutils.exportSettings(settings, path=EXPORT_PATH)
            
            if EXPORT_OUTPUT:
                ioutils.exportOutputs(model_col, dist_col, best_dist, avg_dist, 
                                      settings, t2-t1, ens_st, path=EXPORT_PATH)

        

