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
import matplotlib.pyplot as plt
import time

def f1(k_list, *args):
    global counts
    global countf
    
    args[0].reset()
    
    for i in range(5):
        args[0].model[realGlobalParameterIds[i]] = k_list[i]
        
    ss = args[0].steadyStateSolver
    ss.allow_approx = False
    ss.allow_presimulation = False
    
    try:
        r.steadyState()
        ccTest = args[0].getUnscaledConcentrationControlCoefficientMatrix()
        dist_i = (np.linalg.norm(realConcCC - ccTest))
    except:
        countf += 1
        dist_i = 10000
        
    counts += 1
    
    return dist_i

def callbackF(X, convergence=0.):
    global counts
    global countf
    print(str(counts) + ", " + str(countf))
    return False

def mutate_and_evaluate(listantStr, listdist, mut_ind):
    global countf
    global counts
    
    eval_dist = np.empty(mut_size)
    eval_model = np.empty(mut_size, dtype='object')
    
    for m in mut_range:
        antimony.loadAntimonyString(listantStr[m])
        module = antimony.getModuleNames()[-1]
        
        rct = np.array(antimony.getReactantNames(module)).tolist()
        prd = np.array(antimony.getProductNames(module)).tolist()
        
        r = te.loada(listantStr[m])
        param_val = r.getGlobalParameterValues()
        
        # TODO: Remove part mutation (done)
        # TODO: multiply rate constants by 0.1 , 0.5 etc. (randomize direction for each rate constants?) (done)
        # TODO: randomly choose between  rate constants and reactions (done)
        # TODO: change initial boundary species in similar manner (maybe)
        # TODO: For top performers, only mutate rate constants
        # TODO: assert unique reactions (Done)
        stt = [[],[],[]]
        st = ens_st[0]
        
        n = 0
        
        while (st in ens_st and n < maxIter_mut):
            while len(stt[1]) != realNumFloating or len(stt[2]) != realNumBoundary:
                r_idx = np.random.randint(0, nr)
                rt = ng.pickReactionType()
                if rt == ng.TReactionType.UNIUNI:
                    # UniUni
                    rct_prd = np.random.choice(np.arange(ns), size=2, replace=False)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == ['S'+str(rct_prd[0])]]
                    all_prd = [i for i, x in enumerate(prd) if x == ['S'+str(rct_prd[1])]]
                    
                    while len(set(all_rct) & set(all_prd)) > 0:
                        rct_prd = np.random.choice(np.arange(ns), size=2, replace=False)
                        all_rct = [i for i, x in enumerate(rct) if x == ['S'+str(rct_prd[0])]]
                        all_prd = [i for i, x in enumerate(prd) if x == ['S'+str(rct_prd[1])]]
                    
                    rct[r_idx] = ["S" + str(rct_prd[0])]
                    prd[r_idx] = ["S" + str(rct_prd[1])]
                    
                elif rt == ng.TReactionType.BIUNI:
                    # BiUni
                    # Pick two reactants
                    rct_id = np.random.choice(np.arange(ns), size=2, replace=True)
                    # pick a product but only products that don't include the reactants
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                    
                     # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == ('S'+str(rct_id[0]),
                                                       'S'+str(rct_id[1])) or x == ('S'+str(rct_id[1]), 'S'+str(rct_id[0]))]
                    all_prd = [i for i, x in enumerate(prd) if x == ['S'+str(prd_id[0])]]
                    
                    while len(set(all_rct) & set(all_prd)) > 0:
                        rct_id = np.random.choice(np.arange(ns), size=2, replace=True)
                        prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                        all_rct = [i for i, x in enumerate(rct) if x == ('S'+str(rct_id[0]),
                                                       'S'+str(rct_id[1])) or x == ('S'+str(rct_id[1]), 'S'+str(rct_id[0]))]
                        all_prd = [i for i, x in enumerate(prd) if x == ['S'+str(prd_id[0])]]
                    
                    rct[r_idx] = ["S" + str(rct_id[0]), "S" + str(rct_id[1])]
                    prd[r_idx] = ["S" + str(prd_id[0])]
                    
                elif rt == ng.TReactionType.UNIBI:
                    # UniBi
                    rct_id = np.random.choice(np.arange(ns), size=1)
                    # pick a product but only products that don't include the reactant
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, replace=True)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == ['S'+str(rct_id[0])]]
                    all_prd = [i for i, x in enumerate(prd) if x == ('S'+str(prd_id[0]),
                                                       'S'+str(prd_id[1])) or x == ('S'+str(prd_id[1]), 'S'+str(prd_id[0]))]
                    
                    while len(set(all_rct) & set(all_prd)) > 0:
                        rct_id = np.random.choice(np.arange(ns), size=1)
                        prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, replace=True)
                        all_rct = [i for i, x in enumerate(rct) if x == ['S'+str(rct_id[0])]]
                        all_prd = [i for i, x in enumerate(prd) if x == ('S'+str(prd_id[0]),
                                                       'S'+str(prd_id[1])) or x == ('S'+str(prd_id[1]), 'S'+str(prd_id[0]))]
                    
                    rct[r_idx] = ["S" + str(rct_id[0])]
                    prd[r_idx] = ["S" + str(prd_id[0]), "S" + str(prd_id[1])]
                    
                else:
                    # BiBi
                    rct_id = np.random.choice(np.arange(ns), size=2, replace=True)
                    # pick a product but only products that don't include the reactant
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, replace=True)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == ('S'+str(rct_id[0]),
                                                       'S'+str(rct_id[1])) or x == ('S'+str(rct_id[1]), 'S'+str(rct_id[0]))]
                    all_prd = [i for i, x in enumerate(prd) if x == ('S'+str(prd_id[0]),
                                                       'S'+str(prd_id[1])) or x == ('S'+str(prd_id[1]), 'S'+str(prd_id[0]))]
                    
                    while len(set(all_rct) & set(all_prd)) > 0:
                        rct_id = np.random.choice(np.arange(ns), size=2, replace=True)
                        prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, replace=True)
                        all_rct = [i for i, x in enumerate(rct) if x == ('S'+str(rct_id[0]),
                                                       'S'+str(rct_id[1])) or x == ('S'+str(rct_id[1]), 'S'+str(rct_id[0]))]
                        all_prd = [i for i, x in enumerate(prd) if x == ('S'+str(prd_id[0]),
                                                       'S'+str(prd_id[1])) or x == ('S'+str(prd_id[1]), 'S'+str(prd_id[0]))]
                    
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
                    reactionList.append ([rtype, rct_ind, prd_ind, param_val[i]]) 
                
                reactionList.insert(0, ns)
                
                st = ng.getFullStoichiometryMatrix(reactionList).tolist()
                stt = ng.removeBoundaryNodes(np.array(st))
            n += 1
        
        if n >= maxIter_mut:
            eval_dist[m] = listdist[m]
            eval_model[m] = listantStr[m]
        else:
            antStr = ng.genAntimonyScript(realFloatingIds, realBoundaryIds, stt[1], stt[2], reactionList, boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)
                if (realFloatingIds == np.sort(r.getFloatingSpeciesIds())).all() and (realBoundaryIds == np.sort(r.getBoundarySpeciesIds())).all():
                    
                    counts = 0
                    countf = 0
                    
                    res = scipy.optimize.differential_evolution(f1, args=(r,), bounds=p_bound, 
                                                polish=False, tol=0.03)
                
                    r.reset()
                    
                    for i in range(5):
                        r.model[realGlobalParameterIds[i]] = res.x[i]
                    
                    ss = r.steadyStateSolver
                    ss.allow_approx = False
                    ss.allow_presimulation = False
                    r.steadyState()
                    SS_i = r.getFloatingSpeciesConcentrations()
                    if np.any(SS_i > 1e5):
                        r.reset()
                        ss.allow_presimulation = True
                        ss.presimulation_time = 1000
                        r.steadyState()
                        SS_i = r.getFloatingSpeciesConcentrations()
                    if np.any(SS_i < 1E-5) or np.any(SS_i > 1e5):
                        antimony.clearPreviousLoads()
    #                    if mut_ind[m] < pass_size:
    #                        rnd_dist, rnd_model = random_gen(1)
    #                        eval_dist[m] = rnd_dist[0]
    #                        eval_model[m] = rnd_model[0]
    #                    else:
                        eval_dist[m] = listdist[m]
                        eval_model[m] = listantStr[m]
                    else:
    #                    F_i = sr.getReactionRates()
                        
    #                    fluxCC_i = r.getScaledFluxControlCoefficientMatrix()
    #                    fluxCC_i_row = fluxCC_i.rownames
    #                    fluxCC_i_col = fluxCC_i.colnames
    #                    fluxCC_i = fluxCC_i[np.argsort(fluxCC_i_row)]
    #                    fluxCC_i = fluxCC_i[:,np.argsort(fluxCC_i_col)]
                        
                        concCC_i = r.getUnscaledConcentrationControlCoefficientMatrix()
                        concCC_i_row = concCC_i.rownames
                        concCC_i_col = concCC_i.colnames
                        concCC_i = concCC_i[np.argsort(concCC_i_row)]
                        concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                        
    #                    dist_i = (w1*(np.linalg.norm(realFluxCC - fluxCC_i) + np.linalg.norm(realConcCC - concCC_i))
    #                            + w2*(np.sqrt(np.sum(np.square(realFlux - F_i))) + np.sqrt(np.sum(np.square(realSteadyState - SS_i)))))
                        dist_i = (w1*(np.linalg.norm(realConcCC - concCC_i)))
                        
                        antimony.clearPreviousLoads()
                        if dist_i < listdist[m]:
                            eval_dist[m] = dist_i
                            eval_model[m] = r.getAntimony(current=True)
                            ens_st.append(st)
                        else:
                            eval_dist[m] = listdist[m]
                            eval_model[m] = listantStr[m]
                else:
                    antimony.clearPreviousLoads()
    #                if mut_ind[m] < pass_size:
    #                    rnd_dist, rnd_model = random_gen(1)
    #                    eval_dist[m] = rnd_dist[0]
    #                    eval_model[m] = rnd_model[0]
    #                else:
                    eval_dist[m] = listdist[m]
                    eval_model[m] = listantStr[m]
            except:
                antimony.clearPreviousLoads()
    #            if mut_ind[m] < pass_size:
    #                rnd_dist, rnd_model = random_gen(1)
    #                eval_dist[m] = rnd_dist[0]
    #                eval_model[m] = rnd_model[0]
    #            else:
                eval_dist[m] = listdist[m]
                eval_model[m] = listantStr[m]

    print("Unique models: " + str(len(set(listantStr))))
    print("Unique models: " + str(len(set(eval_model))))

    return eval_dist, eval_model

# TODO: assert the same number of bounary and floating species (Done)
# TODO: Weighting for control coefficients and steady-states/fluxes? (Done)
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
    while numGoodModels < ens_size:
        #antStr = ng.generateLinearChainAnt(ns)
        rl = ng.generateReactionList(ns, nr)
        st = ng.getFullStoichiometryMatrix(rl).tolist()
        # Ensure no redundant model
        while st in ens_st:
            rl = ng.generateReactionList(ns, nr)
            st = ng.getFullStoichiometryMatrix(rl).tolist()
        stt = ng.removeBoundaryNodes(np.array(st))
        if len(stt[1]) == realNumFloating and len(stt[2]) == realNumBoundary:
            antStr = ng.genAntimonyScript(realFloatingIds, realBoundaryIds, stt[1], stt[2], rl, boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)

                counts = 0
                countf = 0
                
                res = scipy.optimize.differential_evolution(f1, args=(r,), bounds=p_bound,
                                                polish=False, tol=0.03)
                r.reset()
                for i in range(5):
                    r.model[realGlobalParameterIds[i]] = res.x[i]
                    
                ss = r.steadyStateSolver
                ss.allow_approx = False
                ss.allow_presimulation = False
                
                r.steadyState()
                SS_i = r.getFloatingSpeciesConcentrations()
                # Buggy model
                if np.any(SS_i > 1e5):
                    r.reset()
                    ss.allow_presimulation = True
                    ss.presimulation_time = 1000
                    r.steadyState()
                    SS_i = r.getFloatingSpeciesConcentrations()
                if np.any(SS_i < 1E-5) or np.any(SS_i > 1e5):
                    continue
                else:
                    #F_i = r.getReactionRates()
                    
#                    fluxCC_i = r.getScaledFluxControlCoefficientMatrix()
#                    fluxCC_i_row = fluxCC_i.rownames
#                    fluxCC_i_col = fluxCC_i.colnames
#                    fluxCC_i = fluxCC_i[np.argsort(fluxCC_i_row)]
#                    fluxCC_i = fluxCC_i[:,np.argsort(fluxCC_i_col)]
                    
                    concCC_i = r.getUnscaledConcentrationControlCoefficientMatrix()
                    concCC_i_row = concCC_i.rownames
                    concCC_i_col = concCC_i.colnames
                    concCC_i = concCC_i[np.argsort(concCC_i_row)]
                    concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                    
#                    dist_i = (w1*(np.linalg.norm(realFluxCC - fluxCC_i) + np.linalg.norm(realConcCC - concCC_i))
#                        + w2*(np.sqrt(np.sum(np.square(realFlux - F_i))) + np.sqrt(np.sum(np.square(realSteadyState - SS_i)))))
                    dist_i = (w1*(np.linalg.norm(realConcCC - concCC_i)))
                    ens_dist[numGoodModels] = dist_i
                    ens_model[numGoodModels] = r.getAntimony(current=True)
                    ens_st.append(st)#[numGoodModels] = st
                    
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


def random_gen(rnd_size):
    global countf
    global counts
    
    numGoodModels = 0
    
    rnd_dist = np.empty(rnd_size)
    rnd_model = np.empty(rnd_size, dtype='object')
    
    while numGoodModels < rnd_size:
        rl = ng.generateReactionList(ns, nr)
        st = ng.getFullStoichiometryMatrix(rl).tolist()
        # Ensure no redundant models
        while st in ens_st:
            rl = ng.generateReactionList(ns, nr)
            st = ng.getFullStoichiometryMatrix(rl).tolist()
        stt = ng.removeBoundaryNodes(np.array(st))
        if len(stt[1]) == realNumFloating and len(stt[2]) == realNumBoundary:
            antStr = ng.genAntimonyScript(realFloatingIds, realBoundaryIds, stt[1], stt[2], rl, boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)
                
                counts = 0
                countf = 0
                
                res = scipy.optimize.differential_evolution(f1, args=(r,), bounds=p_bound,
                                                polish=False, tol=0.03)
                
                r.reset()
                
                for i in range(5):
                    r.model[realGlobalParameterIds[i]] = res.x[i]
                
                ss = r.steadyStateSolver
                ss.allow_approx = False
                ss.allow_presimulation = False
                r.steadyState()
                SS_i = r.getFloatingSpeciesConcentrations()
                if np.any(SS_i > 1e5):
                    r.reset()
                    ss.allow_presimulation = True
                    ss.presimulation_time = 1000
                    r.steadyState()
                    SS_i = r.getFloatingSpeciesConcentrations()
                if np.any(SS_i < 1E-5) or np.any(SS_i > 1e5):
                    continue
                else:
#                    F_i = r.getReactionRates()
                    
#                    fluxCC_i = r.getScaledFluxControlCoefficientMatrix()
#                    fluxCC_i_row = fluxCC_i.rownames
#                    fluxCC_i_col = fluxCC_i.colnames
#                    fluxCC_i = fluxCC_i[np.argsort(fluxCC_i_row)]
#                    fluxCC_i = fluxCC_i[:,np.argsort(fluxCC_i_col)]
                    
                    concCC_i = r.getUnscaledConcentrationControlCoefficientMatrix()
                    concCC_i_row = concCC_i.rownames
                    concCC_i_col = concCC_i.colnames
                    concCC_i = concCC_i[np.argsort(concCC_i_row)]
                    concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                    
#                    dist_i = (w1*(np.linalg.norm(realFluxCC - fluxCC_i) + np.linalg.norm(realConcCC - concCC_i))
#                        + w2*(np.sqrt(np.sum(np.square(realFlux - F_i))) + np.sqrt(np.sum(np.square(realSteadyState - SS_i)))))
                    dist_i = (w1*(np.linalg.norm(realConcCC - concCC_i)))
                    
                    rnd_dist[numGoodModels] = dist_i
                    rnd_model[numGoodModels] = r.getAntimony(current=True)
                    ens_st.append(st)
                    
                    numGoodModels = numGoodModels + 1
            except:
                continue
            antimony.clearPreviousLoads()
    
    return rnd_dist, rnd_model

# TODO: no random network (Done)
# TODO: Record seed, etc. (later on)
# TODO: Exclude best performing model from mutation (Done)
# TODO: change selection process (Pick pair and choose better one) (Done)
# TODO: simulated annealing (multiply to fitness for rate constants)
# TODO: Identical distance is possible even though the model is different. Add bindary counting to distance?
if __name__ == '__main__':
    roadrunner.Config.setValue(roadrunner.Config.ROADRUNNER_DISABLE_WARNINGS, 3)

#%% Settings
    model_type = 'FFL' # 'FFL', 'Linear', 'Nested', 'Branched'
    
    n_gen = 4 # Number of generations
    ens_size = 100 # Size of output ensemble
    pass_size = int(ens_size/10) #20 # Size of models passed on the next generation without mutation
    mut_size = int(ens_size/2) #100
    nrnd_size = pass_size+mut_size
    rnd_size = ens_size-pass_size-mut_size
    
    maxIter_gen = 1000
    maxIter_mut = 1000
    
#    MKP = 0. # Probability of changing rate constants on mutation. Probability of changing reaction is 1 - MKP
#    rateStep = 0.1 # Stepsize for mutating rate constants. Actual stepsize is rateConstant*rateStep
    w1 = 1.5 # Weight for control coefficients when calculating the distance
    w2 = 1.0 # Weight for steady-state and flux when calculating the distance
    
    r_seed = 123123 # random see
    r_roulette = False # Flag for using random roulette or best of pair for selection process
    NOISE = False # Flag for adding Gaussian noise to steady-state and control coefficiant values
    noise_std = 0.1 # Standard deviation of Gaussian noise
    
    PLOT = False # Flag for plots
    SAVE = False # Flag for saving plots

#%%
    if model_type == 'Linear':
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
    elif model_type == 'Nested':
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
    elif model_type == 'FFL':
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
    elif model_type == 'Branched':
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
    
    realRR = te.loada(realModel)
    
    realSS = realRR.steadyStateSolver
    realSS.allow_approx = False
    realSS.allow_presimulation = False
    realSS.presimulation_time = 10000
    
    realNumBoundary = realRR.getNumBoundarySpecies()
    realNumFloating = realRR.getNumFloatingSpecies()
    realFloatingIds = np.sort(realRR.getFloatingSpeciesIds())
    realBoundaryIds = np.sort(realRR.getBoundarySpeciesIds())
    realBoundaryVal = realRR.getBoundarySpeciesConcentrations()
    realGlobalParameterIds = realRR.getGlobalParameterIds()
    
    realRR.steadyState()
    realSteadyState = realRR.getFloatingSpeciesConcentrations()
    realFlux = realRR.getReactionRates()
    realRR.reset()
    realFluxCC = realRR.getScaledFluxControlCoefficientMatrix()
    realConcCC = realRR.getUnscaledConcentrationControlCoefficientMatrix()
    
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
    
#%%
    np.random.seed(r_seed)
    
    best_dist = []
    avg_dist = []
    
    n_range = range(1, n_gen)
    ens_range = range(ens_size)
    mut_range = range(mut_size)
    r_range = range(nr)
    
    p_bound = [(1e-3, 1.)]*nr
    
#%%
    t1 = time.time()
    
    ens_dist, ens_model, ens_st = initialize()
    
    dist_top_ind = np.argsort(ens_dist)
    dist_top = ens_dist[dist_top_ind]
    model_top = ens_model[dist_top_ind]
    
    print("Minimum distance: " + str(dist_top[0]))
    print("Average distance: " + str(np.average(dist_top)))
    print("Unique models: " + str(len(set(ens_model))))
    best_dist.append(dist_top[0])
    avg_dist.append(np.average(dist_top))
    
    # TODO: Remove for loop
    for n in n_range:
        if r_roulette:
            ens_inv = np.divide(1, dist_top[pass_size:])
            ens_prob = np.divide(ens_inv, np.sum(ens_inv))
            mut_ind = np.random.choice(np.arange(pass_size, ens_size), size=mut_size, replace=False, p=ens_prob)
        else:
            mut_p = 1/ens_dist/np.sum(1/ens_dist)
            if mut_size > (ens_size)/2.:
                mut_pair = np.random.choice(np.arange(ens_size), size=(mut_size, 2), replace=True)
            else:
                mut_ind = np.random.choice(np.arange(ens_size), size=(mut_size, 1), replace=False, p=mut_p)
           # mut_ind = np.sort(mut_pair[np.arange(mut_size), np.argmin(dist_top[mut_pair], axis=1)])

#        ens_model[:pass_size] = model_top[:pass_size]
#        ens_dist[:pass_size] = dist_top[:pass_size]
        
        eval_dist, eval_model = mutate_and_evaluate(ens_model[mut_ind], ens_dist[mut_ind], mut_ind)
        
        ens_model[mut_ind] = eval_model
        ens_dist[mut_ind] = eval_dist
        
        print("Unique models post mutation: " + str(len(set(ens_model))))
        
        rnd_dist, rnd_model = random_gen(rnd_size)
        ens_model[nrnd_size:] = rnd_model
        ens_dist[nrnd_size:] = rnd_dist
        
        #print("Unique models post random gen: " + str(len(set(ens_model))))
        
        dist_top_ind = np.argsort(ens_dist)
        dist_top = ens_dist[dist_top_ind]
        model_top = ens_model[dist_top_ind]
        
        print("In generation: " + str(n+1))
        print("Minimum distance: " + str(dist_top[0]))
        print("Average distance: " + str(np.average(dist_top)))
        best_dist.append(dist_top[0])
        avg_dist.append(np.average(dist_top))
        
        # Error check
        #if np.average(dist_top) > 10000:
            #break
                
    print(time.time() - t1)
        
#%%
    if PLOT:
        # Convergence
        plt.plot(best_dist)
        #plt.plot(avg_dist)
        plt.xlabel("Generations", fontsize=15)
        plt.ylabel("Distance", fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        if SAVE:
            plt.savefig(os.path.join('./convergence_' + model_type + '.pdf'), bbox_inches='tight')
        plt.show()
        # TODO: Add polishing with fast optimizer 
            
        # Average residual
        r_real = te.loada(realModel)
        result_real = r_real.simulate(0, 100, 100)
        
        top_result = []
        top_diff = []
        
        for i in ens_range:
            r = te.loada(ens_model[np.argsort(ens_dist)[i]])
            top_sim = r.simulate(0, 100, 100)
            top_result.append(top_sim)
            top_diff.append(np.subtract(result_real[:,1:], top_sim[:,1:]))
    
        percentage = 0.1#float(pass_size)/ens_size
        
        ave_diff = np.average(top_diff[:int(ens_size*percentage)], axis=0)
        
        plt.plot(ave_diff)
        plt.xlabel("Time (s)", fontsize=15)
        plt.ylabel("Residual", fontsize=15)
        plt.legend(["S1","S2","S3"])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        if SAVE:
            plt.savefig(os.path.join('./average_residual_' + model_type + '.pdf'), bbox_inches='tight')
        plt.show()
        
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
                top_diff_k.append(np.sqrt(np.divide(np.sum(np.square(np.subtract(k_real, top_k))),len(k_real))))
            except:
                top_diff_k.append(np.sqrt(np.divide(np.sum(np.square(np.subtract(k_real, top_k[1:]))),len(k_real))))
        
        krmse = top_diff_k[:int(ens_size*percentage)]
        
        plt.hist(krmse, bins=15, density=True)
        plt.xlabel("RMSE", fontsize=15)
        plt.ylabel("Normalized Frequency", fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        if SAVE:
            plt.savefig(os.path.join('./parameter_rmse_' + model_type + '.pdf'), bbox_inches='tight')
        plt.show()

#%%


