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
import networkGenerator as ng
import matplotlib.pyplot as plt
import time

def mutate_and_evaluate(listantStr, listdist, listst, mut_ind):
    eval_dist = np.empty(mut_size)
    eval_model = np.empty(mut_size, dtype='object')
    eval_st = np.empty(mut_size, dtype='object')
    
    for m in mut_range:
        antimony.loadAntimonyString(listantStr[m])
        module = antimony.getModuleNames()[-1]
        
        rct = np.array(antimony.getReactantNames(module)).tolist()
        prd = np.array(antimony.getProductNames(module)).tolist()
        
        r = te.loada(listantStr[m])
        param_id = r.getGlobalParameterIds()
        param_val = r.getGlobalParameterValues()
        flt_id = r.getFloatingSpeciesIds()
        bnd_id = r.getBoundarySpeciesIds()
        bnd_val = r.getBoundarySpeciesConcentrations()
        spe_id = sorted(flt_id + bnd_id)
        
        # TODO: Remove part mutation (done)
        # TODO: multiply rate constants by 0.1 , 0.5 etc. (randomize direction for each rate constants?) (done)
        # TODO: randomly choose between  rate constants and reactions (done)
        # TODO: change initial boundary species in similar manner (maybe)
        # TOSO: For top performers, only mutate rate constants
        if mut_ind[m] < pass_size:
            chprob = 0
        else:
            chprob = MKP
        
        w_mut = np.random.random()
        
        if w_mut < chprob: # Change rate constants
            k_idx = np.random.randint(0, len(param_id))
            if np.random.random() < 0.5: # Increase
                param_val[k_idx] = param_val[k_idx] + param_val[k_idx]*rateStep
                if param_val[k_idx] > 1.:
                    param_val[k_idx] = 1.0
            else: # Decrease
                param_val[k_idx] = param_val[k_idx] - param_val[k_idx]*rateStep
                if param_val[k_idx] < 1e-3:
                    param_val[k_idx] = 1e-3
                    
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
            
            flt_id = [s.replace('S', '') for s in flt_id]
            bnd_id = [s.replace('S', '') for s in bnd_id]
            
            antStr = ng.genAntimonyScript(flt_id, bnd_id, flt_id, bnd_id, reactionList, boundary_init=realBoundaryVal)
            
        # TODO: assert unique reactions (Done)
        elif w_mut >= chprob: # Change reactions
            stt = [[],[],[]]
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
                
                st = ng.getFullStoichiometryMatrix(reactionList)
                stt = ng.removeBoundaryNodes(st)
                
            antStr = ng.genAntimonyScript(realFloatingIds, realBoundaryIds, stt[1], stt[2], reactionList, boundary_init=realBoundaryVal)
        
        try:
            r = te.loada(antStr)
            if (realFloatingIds == np.sort(r.getFloatingSpeciesIds())).all() and (realBoundaryIds == np.sort(r.getBoundarySpeciesIds())).all():
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
                    eval_st[m] = listst[m]
                else:
                    F_i = r.getReactionRates()
                    
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
                    
                    eval_dist[m] = dist_i
                    eval_model[m] = antStr
            else:
                antimony.clearPreviousLoads()
#                if mut_ind[m] < pass_size:
#                    rnd_dist, rnd_model = random_gen(1)
#                    eval_dist[m] = rnd_dist[0]
#                    eval_model[m] = rnd_model[0]
#                else:
                eval_dist[m] = listdist[m]
                eval_model[m] = listantStr[m]
                eval_st[m] = listst[m]
        except:
            antimony.clearPreviousLoads()
#            if mut_ind[m] < pass_size:
#                rnd_dist, rnd_model = random_gen(1)
#                eval_dist[m] = rnd_dist[0]
#                eval_model[m] = rnd_model[0]
#            else:
            eval_dist[m] = listdist[m]
            eval_model[m] = listantStr[m]
            eval_st[m] = listst[m]

    return eval_dist, eval_model, eval_st


def polish(listantStr, listdist):
    polish_dist = np.empty(polish_size-1)
    polish_model = np.empty(polish_size-1, dtype='object')
    
    for m in psize_range:
        antimony.loadAntimonyString(listantStr[m])
        module = antimony.getModuleNames()[-1]
        
        rct = np.array(antimony.getReactantNames(module)).tolist()
        prd = np.array(antimony.getProductNames(module)).tolist()
        
        r = te.loada(listantStr[m])
        param_id = r.getGlobalParameterIds()
        param_val = r.getGlobalParameterValues()
        flt_id = r.getFloatingSpeciesIds()
        bnd_id = r.getBoundarySpeciesIds()
        bnd_val = r.getBoundarySpeciesConcentrations()
        spe_id = sorted(flt_id + bnd_id)
        
        k_idx = np.random.randint(0, len(param_id))
        if np.random.random() < 0.5:
            param_val[k_idx] = param_val[k_idx] + param_val[k_idx]*rateStep
            if param_val[k_idx] > 1.:
                param_val[k_idx] = 1.0
        else:
            param_val[k_idx] = param_val[k_idx] - param_val[k_idx]*rateStep
            if param_val[k_idx] < 1e-3:
                param_val[k_idx] = 1e-3
                
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
        
        reactionList.insert (0, ns)
        
        flt_id_s = [s.replace('S', '') for s in flt_id]
        bnd_id_s = [s.replace('S', '') for s in bnd_id]
        
        antStr = ng.genAntimonyScript(flt_id, bnd_id, flt_id_s, bnd_id_s, reactionList, boundary_init=realBoundaryVal)

        try:
            r = te.loada(antStr)
            if (realFloatingIds == np.sort(r.getFloatingSpeciesIds())).all() and (realBoundaryIds == np.sort(r.getBoundarySpeciesIds())).all():
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
                    polish_dist[m] = listdist[m]
                    polish_model[m] = listantStr[m]
                else:
                    F_i = r.getReactionRates()
                    dist_i = w2*(np.sqrt(np.sum(np.square(realFlux - F_i))) + np.sqrt(np.sum(np.square(realSteadyState - SS_i))))
                    #np.sqrt(np.sum(np.square(realSteadyState - SS_i)))
                    
                    antimony.clearPreviousLoads()
                    
                    polish_dist[m] = dist_i
                    polish_model[m] = antStr
            else:
                antimony.clearPreviousLoads()
                polish_dist[m] = listdist[m]
                polish_model[m] = listantStr[m]
        except:
            antimony.clearPreviousLoads()
            polish_dist[m] = listdist[m]
            polish_model[m] = listantStr[m]

    return polish_dist, polish_model

    
# TODO: assert the same number of bounary and floating species (Done)
# TODO: Weighting for control coefficients and steady-states/fluxes? (Done)
def initialize():
    numBadModels = 0
    numGoodModels = 0
    numIter = 0
    
    ens_dist = np.empty(ens_size)
    ens_model = np.empty(ens_size, dtype='object')
    ens_st = np.empty(ens_size, dtype='object')
    
    # Initial Random generation
    while numGoodModels < ens_size:
        #antStr = ng.generateLinearChainAnt(ns)
        rl = ng.generateReactionList(ns, nr)
        st = ng.getFullStoichiometryMatrix(rl)
        # Ensure no redundant model
        while st in ens_st.any():
            rl = ng.generateReactionList(ns, nr)
            st = ng.getFullStoichiometryMatrix(rl)
        stt = ng.removeBoundaryNodes(st)
        if len(stt[1]) == realNumFloating and len(stt[2]) == realNumBoundary:
            antStr = ng.genAntimonyScript(realFloatingIds, realBoundaryIds, stt[1], stt[2], rl, boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)
                #if (realFloatingIds == np.sort(r.getFloatingSpeciesIds())).all() and (realBoundaryIds == np.sort(r.getBoundarySpeciesIds())).all():
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
                    ens_model[numGoodModels] = antStr
                    ens_st[numGoodModels] = st
                    
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
    numGoodModels = 0
    
    rnd_dist = np.empty(rnd_size)
    rnd_model = np.empty(rnd_size, dtype='object')
    
    while numGoodModels < rnd_size:
        rl = ng.generateReactionList(ns, nr)
        st = ng.getFullStoichiometryMatrix(rl)
        stt = ng.removeBoundaryNodes(st)
        if len(stt[1]) == realNumFloating and len(stt[2]) == realNumBoundary:
            antStr = ng.genAntimonyScript(realFloatingIds, realBoundaryIds, stt[1], stt[2], rl, boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)
                #if (realFloatingIds == np.sort(r.getFloatingSpeciesIds())).all() and (realBoundaryIds == np.sort(r.getBoundarySpeciesIds())).all():
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
                    F_i = r.getReactionRates()
                    
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
                    rnd_model[numGoodModels] = antStr
                    
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
    
    n_gen = 100 # Number of generations
    ens_size = 100 # Size of output ensemble
    pass_size = int(ens_size/10) #20 # Size of models passed on the next generation without mutation
    mut_size = int(ens_size/2) #100
    n_polish = 50 # Number of polishing steps to run
    polish_size = 100 # Size of polish population
    nrnd_size = pass_size+mut_size
    rnd_size = ens_size-pass_size-mut_size
    
    MKP = 0.3 # Probability of changing rate constants on mutation. Probability of changing reaction is 1 - MKP
    rateStep = 0.1 # Stepsize for mutating rate constants. Actual stepsize is rateConstant*rateStep
    w1 = 1.5 # Weight for control coefficients when calculating the distance
    w2 = 1.0 # Weight for steady-state and flux when calculating the distance
    
    r_seed = 123123 # random see
    r_roulette = False # Flag for using random roulette or best of pair for selection process
    NOISE = False # Flag for adding Gaussian noise to steady-state and control coefficiant values
    noise_std = 0.1 # Standard deviation of Gaussian noise
    
    PLOT = True # Flag for plots
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
    p_range = range(0, n_polish)
    ens_range = range(ens_size)
    mut_range = range(mut_size)
    psize_range = range(polish_size-1)
    r_range = range(nr)
    
#%%
    t1 = time.time()
    
    ens_dist, ens_model, ens_st = initialize()
    
    dist_top_ind = np.argsort(ens_dist)
    dist_top = ens_dist[dist_top_ind]
    model_top = ens_model[dist_top_ind]
    st_top = ens_st[dist_top_ind]
    
    print("Minimum distance: " + str(dist_top[0]))
    print("Average distance: " + str(np.average(dist_top)))
    best_dist.append(dist_top[0])
    avg_dist.append(np.average(dist_top))
    
    # TODO: Remove for loop
    for n in n_range:
        if r_roulette:
            ens_inv = np.divide(1, dist_top[pass_size:])
            ens_prob = np.divide(ens_inv, np.sum(ens_inv))
            mut_ind = np.random.choice(np.arange(pass_size, ens_size), size=mut_size, replace=False, p=ens_prob)
        else:
            if mut_size > (ens_size)/2.:
                mut_pair = np.random.choice(np.arange(ens_size), size=(mut_size, 2), replace=True)
            else:
                mut_pair = np.random.choice(np.arange(ens_size), size=(mut_size, 2), replace=False)
            mut_ind = np.sort(mut_pair[np.arange(mut_size), np.argmin(dist_top[mut_pair], axis=1)])

        ens_model[:pass_size] = model_top[:pass_size]
        ens_dist[:pass_size] = dist_top[:pass_size]
        ens_st[:pass_size] = st_top[:pass_size]
        
        eval_dist, eval_model, eval_st = mutate_and_evaluate(model_top[mut_ind], dist_top[mut_ind], st_top[mut_ind], mut_ind)
        
        ens_model[pass_size:nrnd_size] = eval_model
        ens_dist[pass_size:nrnd_size] = eval_dist
        ens_st[pass_size:nrnd_size] = eval_st
        
        rnd_dist, rnd_model, rnd_st = random_gen(rnd_size)
        ens_model[nrnd_size:] = rnd_model
        ens_dist[nrnd_size:] = rnd_dist
        ens_st[nrnd_size:] = rnd_st
        
        dist_top_ind = np.argsort(ens_dist)
        dist_top = ens_dist[dist_top_ind]
        model_top = ens_model[dist_top_ind]
        st_top = ens_st[dist_top_ind]
        
        print("In generation: " + str(n+1))
        print("Minimum distance: " + str(dist_top[0]))
        print("Average distance: " + str(np.average(dist_top)))
        best_dist.append(dist_top[0])
        avg_dist.append(np.average(dist_top))
        
        # Error check
        #if np.average(dist_top) > 10000:
            #break
                
    print(time.time() - t1)
    
#%% Polishing
    
    pol_dist = np.empty(polish_size)
    pol_model = np.empty(polish_size, dtype='object')
    best_pdist = []
    
    r = te.loada(model_top[0])
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
    pdist_init = np.sqrt(np.sum(np.square(realSteadyState - SS_i)))
        
    pmodel_top = np.repeat(model_top[0], polish_size) # Now choose the top model
    pdist_top = np.repeat(pdist_init, polish_size)
    
    for p in p_range:
        mut_pair = np.random.choice(np.arange(polish_size), size=(polish_size-1, 2), replace=True)
        mut_ind = np.sort(mut_pair[np.arange(polish_size-1), np.argmin(pdist_top[mut_pair], axis=1)])
            
        pol_model[0] = pmodel_top[0]
        pol_dist[0] = pdist_top[0]
        
        polish_dist, polish_model = polish(pmodel_top[mut_ind], pdist_top[mut_ind])
        
        pol_model[1:] = polish_model
        pol_dist[1:] = polish_dist
        
        pdist_top = pol_dist[np.argsort(pol_dist)]
        pmodel_top = pol_model[np.argsort(pol_dist)]

        print("In polishing step: " + str(p+1))
        print("Minimum distance: " + str(pdist_top[0]))
        best_pdist.append(pdist_top[0])
    
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
        
        # Polishing convergence
        plt.plot(best_pdist)
        #plt.plot(avg_dist)
        plt.xlabel("Generations", fontsize=15)
        plt.ylabel("Distance", fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        if SAVE:
            plt.savefig(os.path.join('./convergence_p_' + model_type + '.pdf'), bbox_inches='tight')
        plt.show()
            
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
        
        # Best polished residual
        r = te.loada(pol_model[0])
        top_psim = r.simulate(0, 100, 100)
        top_pdiff = np.subtract(result_real[:,1:], top_psim[:,1:])
    
        plt.plot(top_pdiff)
        plt.xlabel("Time (s)", fontsize=15)
        plt.ylabel("Residual", fontsize=15)
        plt.legend(["S1","S2","S3"])
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        if SAVE:
            plt.savefig(os.path.join('./average_residual_p_' + model_type + '.pdf'), bbox_inches='tight')
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


