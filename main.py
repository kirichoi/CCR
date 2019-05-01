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
import copy

#np.seterr(all='raise')

def f1(k_list, *args):
    global counts
    global countf
    
    args[0].reset()
    
    args[0].setValues(args[0].getGlobalParameterIds(), k_list)
    
    try:
        args[0].steadyStateApproximate()
        objCCC = args[0].getScaledConcentrationControlCoefficientMatrix()
        if np.isnan(objCCC).any():
            dist_obj = 10000
        else:
            objCCC[np.abs(objCCC) < 1e-12] = 0 # Set small values to zero
            
            objCCC_row = objCCC.rownames
            objCCC_col = objCCC.colnames
            objCCC = objCCC[np.argsort(objCCC_row)]
            objCCC = objCCC[:,np.argsort(objCCC_col)]
            
    #        testCount = np.array(np.unravel_index(np.argsort(objCC, axis=None), objCC.shape)).T
            dist_obj = w1*((np.linalg.norm(realConcCC - objCCC))/
                        (1 + np.sum(np.equal(np.sign(np.array(realConcCC)), np.sign(np.array(objCCC))))))
                    
    #        + np.sum(args[0].getReactionRates() < 0))
    #                    + 1/(1 + np.sum((testCount == realCount).all(axis=1))))
    #        dist_obj = (np.linalg.norm(realConcCC - objCC))
        
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


def mutate_and_evaluate(listantStr, listdist, listrl):
    global countf
    global counts
    
    eval_dist = np.empty(mut_size)
    eval_model = np.empty(mut_size, dtype='object')
    eval_rl = np.empty(mut_size, dtype='object')
    
    for m in mut_range:
        antimony.loadAntimonyString(listantStr[m])
        module = antimony.getModuleNames()[-1]
        
        r = te.loada(listantStr[m])
        
        r.steadyStateApproximate()
        
        tempdiff = np.max(np.abs(realConcCC - 
                   r.getScaledConcentrationControlCoefficientMatrix()), axis=0)
    
        stt = [[],[],[]]
        reactionList = rl_track[0]
        
        o = 0
        
#        connected = False
    
        while ((stt[1] != realFloatingIdsInd or stt[2] != realBoundaryIdsInd or
                reactionList in rl_track) and (o < maxIter_mut)):
            rct = []
            rcttuple = antimony.getReactantNames(module)
            for rcti in range(len(rcttuple)):
                rct_temp = []
                for rctj in range(len(rcttuple[rcti])):
                    rct_temp.append(int(rcttuple[rcti][rctj][1:]))
                rct.append(rct_temp)
            
            prd = []
            prdtuple = antimony.getProductNames(module)
            for prdi in range(len(prdtuple)):
                prd_temp = []
                for prdj in range(len(prdtuple[prdi])):
                    prd_temp.append(int(prdtuple[prdi][prdj][1:]))
                prd.append(prd_temp)
            
            r_idx = np.random.choice(np.arange(nr), 
                                     p=np.divide(tempdiff,np.sum(tempdiff)))
            posRctInd = np.append(np.array(realFloatingIdsInd)[np.where(
                        np.abs(realConcCC[:,r_idx]) > 1e-12)[0]], np.array(realBoundaryIdsInd)).astype(int)
            rType, regType, revType = ng.pickReactionType()
            reactionList = copy.deepcopy(listrl[m])
            
            if rType == ng.ReactionType.UNIUNI:
                # TODO: pick reactants and products based on control coefficients
                # UniUni
                rct_id = np.random.choice(posRctInd, size=1)
                prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                
                while (((np.any(np.isin(rct_id, realBoundaryIdsInd))) and 
                       (np.any(np.isin(prd_id, realBoundaryIdsInd)))) or
                       (len(set(all_rct) & set(all_prd)) > 0)):
                    rct_id = np.random.choice(posRctInd, size=1)
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                    all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                    
                if regType == ng.RegulationType.DEFAULT:
                    act_id = []
                    inhib_id = []
                elif regType == ng.RegulationType.INHIBITION:
                    act_id = []
                    inhib_id = np.random.choice(np.delete(np.arange(ns), 
                               np.unique(np.concatenate([rct_id, prd_id]))), size=1).tolist()
                    if len(inhib_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                elif regType == ng.RegulationType.ACTIVATION:
                    act_id = np.random.choice(np.delete(np.arange(ns), 
                             np.unique(np.concatenate([rct_id, prd_id]))), size=1).tolist()
                    inhib_id = []
                    if len(act_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                else:
                    reg_id = np.random.choice(np.delete(np.arange(ns), 
                             np.unique(np.concatenate([rct_id, prd_id]))), size=2).tolist()
                    act_id = [reg_id[0]]
                    inhib_id = [reg_id[1]]
                    if len(reg_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                        
                reactionList[r_idx] = [rType, regType, revType, 
                                       [rct_id[0]], [prd_id[0]], 
                                       act_id, inhib_id]
                
            elif rType == ng.ReactionType.BIUNI:
                # BiUni
                rct_id = np.random.choice(posRctInd, size=2, replace=True)
                prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                
                while (((np.any(np.isin(rct_id, realBoundaryIdsInd))) and 
                       (np.any(np.isin(prd_id, realBoundaryIdsInd)))) or
                       (len(set(all_rct) & set(all_prd)) > 0)):
                    rct_id = np.random.choice(posRctInd, size=2, replace=True)
                    # pick a product but only products that don't include the reactants
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=1)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                    all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                    
                if regType == ng.RegulationType.DEFAULT:
                    act_id = []
                    inhib_id = []
                elif regType == ng.RegulationType.INHIBITION:
                    act_id = []
                    inhib_id = np.random.choice(np.delete(np.arange(ns), 
                               np.unique(np.concatenate([rct_id, prd_id]))), size=1).tolist()
                    if len(inhib_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                elif regType == ng.RegulationType.ACTIVATION:
                    act_id = np.random.choice(np.delete(np.arange(ns), 
                             np.unique(np.concatenate([rct_id, prd_id]))), size=1).tolist()
                    inhib_id = []
                    if len(act_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                else:
                    reg_id = np.random.choice(np.delete(np.arange(ns), 
                             np.unique(np.concatenate([rct_id, prd_id]))), size=2).tolist()
                    act_id = [reg_id[0]]
                    inhib_id = [reg_id[1]]
                    if len(reg_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                
                reactionList[r_idx] = [rType, regType, revType,
                                       [rct_id[0], rct_id[1]], [prd_id[0]], 
                                       act_id, inhib_id]
                
                
            elif rType == ng.ReactionType.UNIBI:
                # UniBi
                rct_id = np.random.choice(posRctInd, size=1)
                prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, replace=True)
                all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                
                while (((np.any(np.isin(rct_id, realBoundaryIdsInd))) and 
                       (np.any(np.isin(prd_id, realBoundaryIdsInd)))) or
                       (len(set(all_rct) & set(all_prd)) > 0)):
                    rct_id = np.random.choice(posRctInd, size=1)
                    # pick a product but only products that don't include the reactant
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, replace=True)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                    all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                    
                if regType == ng.RegulationType.DEFAULT:
                    act_id = []
                    inhib_id = []
                elif regType == ng.RegulationType.INHIBITION:
                    act_id = []
                    inhib_id = np.random.choice(np.delete(np.arange(ns), 
                               np.unique(np.concatenate([rct_id, prd_id]))), size=1).tolist()
                    if len(inhib_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                elif regType == ng.RegulationType.ACTIVATION:
                    act_id = np.random.choice(np.delete(np.arange(ns), 
                             np.unique(np.concatenate([rct_id, prd_id]))), size=1).tolist()
                    inhib_id = []
                    if len(act_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                else:
                    reg_id = np.random.choice(np.delete(np.arange(ns), 
                             np.unique(np.concatenate([rct_id, prd_id]))), size=2).tolist()
                    act_id = [reg_id[0]]
                    inhib_id = [reg_id[1]]
                    if len(reg_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                        
                reactionList[r_idx] = [rType, regType, revType,
                                       [rct_id[0]], [prd_id[0], prd_id[1]], 
                                       act_id, inhib_id]
                
            else:
                # BiBi
                rct_id = np.random.choice(posRctInd, size=2, replace=True)
                prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, replace=True)
                all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                
                while (((np.any(np.isin(rct_id, realBoundaryIdsInd))) and 
                       (np.any(np.isin(prd_id, realBoundaryIdsInd)))) or
                       (len(set(all_rct) & set(all_prd)) > 0)):
                    rct_id = np.random.choice(posRctInd, size=2, replace=True)
                    # pick a product but only products that don't include the reactant
                    prd_id = np.random.choice(np.delete(np.arange(ns), rct_id), size=2, replace=True)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                    all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                    
                if regType == ng.RegulationType.DEFAULT:
                    act_id = []
                    inhib_id = []
                elif regType == ng.RegulationType.INHIBITION:
                    act_id = []
                    inhib_id = np.random.choice(np.delete(np.arange(ns), 
                               np.unique(np.concatenate([rct_id, prd_id]))), size=1).tolist()
                    if len(inhib_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                elif regType == ng.RegulationType.ACTIVATION:
                    act_id = np.random.choice(np.delete(np.arange(ns), 
                             np.unique(np.concatenate([rct_id, prd_id]))), size=1).tolist()
                    inhib_id = []
                    if len(act_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                else:
                    reg_id = np.random.choice(np.delete(np.arange(ns), 
                             np.unique(np.concatenate([rct_id, prd_id]))), size=2).tolist()
                    act_id = [reg_id[0]]
                    inhib_id = [reg_id[1]]
                    if len(reg_id) == 0:
                        regType = ng.RegulationType.DEFAULT
                
                reactionList[r_idx] = [rType, regType, revType,
                                       [rct_id[0], rct_id[1]], [prd_id[0], prd_id[1]], 
                                       act_id, inhib_id]
                
#            connected = analysis.isConnected(reactionList)
            st = ng.getFullStoichiometryMatrix(reactionList, ns).tolist()
            stt = ng.removeBoundaryNodes(np.array(st))
            o += 1
        
        if o >= maxIter_mut:
            eval_dist[m] = listdist[m]
            eval_model[m] = listantStr[m]
            eval_rl[m] = listrl[m]
        else:
            antStr = ng.generateAntimony(realFloatingIds, realBoundaryIds, 
                                          stt[1], stt[2], reactionList, 
                                          boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)
                
                counts = 0
                countf = 0
                
                r.steadyStateApproximate()

                p_bound = ng.generateParameterBoundary(r.getGlobalParameterIds())
                res = scipy.optimize.differential_evolution(f1, args=(r,), 
                            bounds=p_bound, maxiter=optiMaxIter, tol=optiTol,
                            polish=optiPolish, seed=r_seed)
                
                if not res.success:
                    eval_dist[m] = listdist[m]
                    eval_model[m] = listantStr[m]
                    eval_rl[m] = listrl[m]
                else:
                    r = te.loada(antStr)
                    r.setValues(r.getGlobalParameterIds(), res.x)
                    
                    r.steadyStateApproximate()
                    SS_i = r.getFloatingSpeciesConcentrations()
                    
                    r.steadyStateApproximate()
                    
                    if np.any(SS_i < 1e-5) or np.any(SS_i > 1e5):
                        eval_dist[m] = listdist[m]
                        eval_model[m] = listantStr[m]
                        eval_rl[m] = listrl[m]
                    else:
                        concCC_i = r.getScaledConcentrationControlCoefficientMatrix()
                        
                        if np.isnan(concCC_i).any():
                            eval_dist[m] = listdist[m]
                            eval_model[m] = listantStr[m]
                            eval_rl[m] = listrl[m]
                        else:
                            concCC_i[np.abs(concCC_i) < 1e-12] = 0 # Set small values to zero
                            
                            concCC_i_row = concCC_i.rownames
                            concCC_i_col = concCC_i.colnames
                            concCC_i = concCC_i[np.argsort(concCC_i_row)]
                            concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                            
#                            count_i = np.array(np.unravel_index(np.argsort(concCC_i, axis=None), concCC_i.shape)).T
                            dist_i = w1*((np.linalg.norm(realConcCC - concCC_i))/
                                      (1 + np.sum(np.equal(np.sign(np.array(realConcCC)), np.sign(np.array(concCC_i))))))
#                            + np.sum(r.getReactionRates() < 0))
#                                        + 1/(1 + np.sum((count_i == realCount).all(axis=1))))
                            
#                            dist_i = w1*(np.linalg.norm(realConcCC - concCC_i))
                            
                            if dist_i < listdist[m]:
                                eval_dist[m] = dist_i
                                r.reset()
                                eval_model[m] = r.getAntimony(current=True)
                                eval_rl[m] = reactionList
                                rl_track.append(reactionList)
                            else:
                                eval_dist[m] = listdist[m]
                                eval_model[m] = listantStr[m]
                                eval_rl[m] = listrl[m]
            except:
                eval_dist[m] = listdist[m]
                eval_model[m] = listantStr[m]
                eval_rl[m] = listrl[m]
        antimony.clearPreviousLoads()

    return eval_dist, eval_model, eval_rl


def initialize():
    global countf
    global counts
    
    numBadModels = 0
    numGoodModels = 0
    numIter = 0
    
    ens_dist = np.empty(ens_size)
    ens_model = np.empty(ens_size, dtype='object')
    ens_rl = np.empty(ens_size, dtype='object')
    rl_track = []
    
    # Initial Random generation
    while (numGoodModels < ens_size):
        rl = ng.generateReactionList(ns, nr, realBoundaryIdsInd)
        st = ng.getFullStoichiometryMatrix(rl, ns).tolist()
        stt = ng.removeBoundaryNodes(np.array(st))
        # Ensure no redundant model
        while (stt[1] != realFloatingIdsInd or stt[2] != realBoundaryIdsInd 
               or rl in rl_track):
            rl = ng.generateReactionList(ns, nr, realBoundaryIdsInd)
            st = ng.getFullStoichiometryMatrix(rl, ns).tolist()
            stt = ng.removeBoundaryNodes(np.array(st))
        antStr = ng.generateAntimony(realFloatingIds, realBoundaryIds, stt[1],
                                      stt[2], rl, boundary_init=realBoundaryVal)
        try:
            r = te.loada(antStr)

            counts = 0
            countf = 0
            
            r.steadyStateApproximate()
            
            p_bound = ng.generateParameterBoundary(r.getGlobalParameterIds())
            res = scipy.optimize.differential_evolution(f1, args=(r,), 
                               bounds=p_bound, maxiter=optiMaxIter, tol=optiTol,
                               polish=optiPolish, seed=r_seed)
            
            if not res.success:
                numBadModels += 1
            else:
                # TODO: Might be able to cut the bottom part by simply using 
                # the obj func value from optimizer
                r = te.loada(antStr)
                r.setValues(r.getGlobalParameterIds(), res.x)
                    
                r.steadyStateApproximate()
                SS_i = r.getFloatingSpeciesConcentrations()
                
                r.steadyStateApproximate()
                
                if np.any(SS_i < 1e-5) or np.any(SS_i > 1e5):
                    numBadModels += 1
                else:
                    concCC_i = r.getScaledConcentrationControlCoefficientMatrix()
                    
                    if np.isnan(concCC_i).any():
                        numBadModels += 1
                    else:
                        concCC_i[np.abs(concCC_i) < 1e-12] = 0 # Set small values to zero
                        
                        concCC_i_row = concCC_i.rownames
                        concCC_i_col = concCC_i.colnames
                        concCC_i = concCC_i[np.argsort(concCC_i_row)]
                        concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                        
#                        count_i = np.array(np.unravel_index(np.argsort(concCC_i, axis=None), concCC_i.shape)).T
                        dist_i = w1*((np.linalg.norm(realConcCC - concCC_i))/
                                      (1 + np.sum(np.equal(np.sign(np.array(realConcCC)), np.sign(np.array(concCC_i))))))
#                        + np.sum(r.getReactionRates() < 0))
#                                    + 1/(1 + np.sum((count_i == realCount).all(axis=1))))
                        
#                        dist_i = w1*(np.linalg.norm(realConcCC - concCC_i))
                        
                        ens_dist[numGoodModels] = dist_i
                        r.reset()
                        ens_model[numGoodModels] = r.getAntimony(current=True)
                        ens_rl[numGoodModels] = rl
                        rl_track.append(rl)
                        
                        print(numGoodModels)
                        print(r.getAntimony(current=True))
                        
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
    
    return ens_dist, ens_model, ens_rl, rl_track


def random_gen(listAntStr, listDist, listrl):
    global countf
    global counts
    
    rndSize = len(listDist)
    
    rnd_dist = np.empty(rndSize)
    rnd_model = np.empty(rndSize, dtype='object')
    rnd_rl = np.empty(rndSize, dtype='object')
    
    for l in range(rndSize):
        d = 0
        rl = ng.generateReactionList(ns, nr, realBoundaryIdsInd)
        st = ng.getFullStoichiometryMatrix(rl, ns).tolist()
        stt = ng.removeBoundaryNodes(np.array(st))
        # Ensure no redundant models
        while ((stt[1] != realFloatingIdsInd or stt[2] != realBoundaryIdsInd or
                rl in rl_track) and (d < maxIter_gen)):
            rl = ng.generateReactionList(ns, nr, realBoundaryIdsInd)
            st = ng.getFullStoichiometryMatrix(rl, ns).tolist()
            stt = ng.removeBoundaryNodes(np.array(st))
            d += 1
            
        if d >= maxIter_gen:
            rnd_dist[l] = listDist[l]
            rnd_model[l] = listAntStr[l]
            rnd_rl[l] = listrl[l]
        else:
            antStr = ng.generateAntimony(realFloatingIds, realBoundaryIds, 
                            stt[1], stt[2], rl, boundary_init=realBoundaryVal)
            try:
                r = te.loada(antStr)
                
                counts = 0
                countf = 0
                
                r.steadyStateApproximate()
                
                p_bound = ng.generateParameterBoundary(r.getGlobalParameterIds())
                res = scipy.optimize.differential_evolution(f1, args=(r,), 
                            bounds=p_bound, maxiter=optiMaxIter, tol=optiTol,
                            polish=optiPolish, seed=r_seed)
                
                # Failed to find solution
                if not res.success:
                    rnd_dist[l] = listDist[l]
                    rnd_model[l] = listAntStr[l]
                    rnd_rl[l] = listrl[l]
                else:
                    r = te.loada(antStr)
                    r.setValues(r.getGlobalParameterIds(), res.x)

                    r.steadyStateApproximate()
                    SS_i = r.getFloatingSpeciesConcentrations()
                    
                    r.steadyStateApproximate()
                    
                    if np.any(SS_i < 1e-5) or np.any(SS_i > 1e5):
                        rnd_dist[l] = listDist[l]
                        rnd_model[l] = listAntStr[l]
                        rnd_rl[l] = listrl[l]
                    else:
                        concCC_i = r.getScaledConcentrationControlCoefficientMatrix()
                        
                        if np.isnan(concCC_i).any():
                            rnd_dist[l] = listDist[l]
                            rnd_model[l] = listAntStr[l]
                            rnd_rl[l] = listrl[l]
                        else:
                            concCC_i[np.abs(concCC_i) < 1e-12] = 0 # Set small values to zero
                            
                            concCC_i_row = concCC_i.rownames
                            concCC_i_col = concCC_i.colnames
                            concCC_i = concCC_i[np.argsort(concCC_i_row)]
                            concCC_i = concCC_i[:,np.argsort(concCC_i_col)]
                            
#                            count_i = np.array(np.unravel_index(np.argsort(concCC_i, axis=None), concCC_i.shape)).T
                            dist_i = w1*((np.linalg.norm(realConcCC - concCC_i))/
                                      (1 + np.sum(np.equal(np.sign(np.array(realConcCC)), np.sign(np.array(concCC_i))))))
#                            + np.sum(r.getReactionRates() < 0))
#                                        + 1/(1 + np.sum((count_i == realCount).all(axis=1))))
                            
#                            dist_i = w1*(np.linalg.norm(realConcCC - concCC_i))
                            
                            if dist_i < listDist[l]:
                                rnd_dist[l] = dist_i
                                r.reset()
                                rnd_model[l] = r.getAntimony(current=True)
                                rnd_rl[l] = rl
                                rl_track.append(rl)
                            else:
                                rnd_dist[l] = listDist[l]
                                rnd_model[l] = listAntStr[l]
                                rnd_rl[l] = listrl[l]
            except:
                rnd_dist[l] = listDist[l]
                rnd_model[l] = listAntStr[l]
                rnd_rl[l] = listrl[l]
        antimony.clearPreviousLoads()
        
    return rnd_dist, rnd_model, rnd_rl

# TODO: simulated annealing (multiply to fitness for rate constants)
if __name__ == '__main__':
    roadrunner.Config.setValue(roadrunner.Config.ROADRUNNER_DISABLE_WARNINGS, 3)


#%% Settings
    
    # Input data
    INPUT = None
    READ_SETTINGS = None
    
    # Test models =============================================================
    
    # 'FFL', 'Linear', 'Nested', 'Branched'
    modelType = 'Full1' 
    
    
    # General settings ========================================================
    
    # Number of generations
    n_gen = 500
    # Size of output ensemble
    ens_size = 200
    # Number of models passed on the next generation without mutation
    pass_size = int(ens_size/10)
    # Number of models to mutate
    mut_size = int(ens_size/2)
    # Maximum iteration allowed for random generation
    maxIter_gen = 20
    # Maximum iteration allowed for mutation
    maxIter_mut = 20
    # Set conserved moiety
    conservedMoiety = True
    
    
    # Optimizer settings ======================================================
    
    # Maximum iteration allowed for optimizer
    optiMaxIter = 100
    optiTol = 1.
    optiPolish = False
    # Weight for control coefficients when calculating the distance
    w1 = 16
    # Weight for steady-state and flux when calculating the distance
    w2 = 1.0
    
    
    # Random settings =========================================================
    
    # random seed
    r_seed = 123123
    # Flag for adding Gaussian noise to steady-state and control coefficiant values
    NOISE = False
    # Standard deviation of Gaussian noise
    ABS_NOISE_STD = 0.01
    # Standard deviation of Gaussian noise
    REL_NOISE_STD = 0.1
    
    
    # Plotting settings =======================================================
    
    # Flag for plots
    PLOT = True
    # Flag for saving plots
    SAVE_PLOT = True
    
    
    # Data settings ===========================================================
    
    # Flag for collecting models
    EXPORT_ALL_MODELS = True
    # Flag for saving collected models
    EXPORT_OUTPUT = True
    # Flag for saving current settings
    EXPORT_SETTINGS = False
    # Path to save the output
    EXPORT_PATH = './USE/output_Full1_3'
    
    # Flag to run algorithm
    RUN = True
    
#%% Analyze True Model
    if conservedMoiety:
        roadrunner.Config.setValue(roadrunner.Config.LOADSBMLOPTIONS_CONSERVED_MOIETIES, True)
    
    roadrunner.Config.setValue(roadrunner.Config.STEADYSTATE_APPROX, True)
    roadrunner.Config.setValue(roadrunner.Config.STEADYSTATE_APPROX_MAX_STEPS, 15)
    roadrunner.Config.setValue(roadrunner.Config.STEADYSTATE_APPROX_TIME, 10000)
#    roadrunner.Config.setValue(roadrunner.Config.STEADYSTATE_APPROX_TOL, 1e-3)
        
    # Using one of the test models
    realModel = ioutils.testModels(modelType)
    
    realRR = te.loada(realModel)
    
    # Species
    realNumBoundary = realRR.getNumBoundarySpecies()
    realNumFloating = realRR.getNumFloatingSpecies()
    realFloatingIds = np.sort(realRR.getFloatingSpeciesIds())
    realFloatingIdsInd = list(map(int, [s.strip('S') for s in realRR.getFloatingSpeciesIds()]))
    realFloatingIdsInd.sort()
    realBoundaryIds = np.sort(realRR.getBoundarySpeciesIds())
    realBoundaryIdsInd = list(map(int,[s.strip('S') for s in realRR.getBoundarySpeciesIds()]))
    realBoundaryIdsInd.sort()
    realBoundaryVal = realRR.getBoundarySpeciesConcentrations()
    realGlobalParameterIds = realRR.getGlobalParameterIds()
    
    # Control Coefficients and Fluxes
    
    realRR.steadyState()
    realSteadyState = realRR.getFloatingSpeciesConcentrations()
    realSteadyStateRatio = np.divide(realSteadyState, np.min(realSteadyState))
    realFlux = realRR.getReactionRates()
    realRR.reset()
    realRR.steadyState()
    realFluxCC = realRR.getScaledFluxControlCoefficientMatrix()
    realConcCC = realRR.getScaledConcentrationControlCoefficientMatrix()
    
    realFluxCC[np.abs(realFluxCC) < 1e-12] = 0
    realConcCC[np.abs(realConcCC) < 1e-12] = 0
    
    # Ordering
    realFluxCCrow = realFluxCC.rownames
    realFluxCCcol = realFluxCC.colnames
    realFluxCC = realFluxCC[np.argsort(realFluxCCrow)]
    realFluxCC = realFluxCC[:,np.argsort(realFluxCCcol)]
    
    realConcCCrow = realConcCC.rownames
    realConcCCcol = realConcCC.colnames
    realConcCC = realConcCC[np.argsort(realConcCCrow)]
    realConcCC = realConcCC[:,np.argsort(realConcCCcol)]
    
    realFlux = realFlux[np.argsort(realRR.getFloatingSpeciesIds())]
    
    # Number of Species and Ranges
    ns = realNumBoundary + realNumFloating # Number of species
    nr = realRR.getNumReactions() # Number of reactions
    
    n_range = range(1, n_gen)
    ens_range = range(ens_size)
    mut_range = range(mut_size)
    r_range = range(nr)
    
    realCount = np.array(np.unravel_index(np.argsort(realFluxCC, axis=None), realFluxCC.shape)).T
        
    #%%
    if RUN:
        
        print("Original Control Coefficients")
        print(realConcCC)
        print("Original Steady State Ratio")
        print(realSteadyStateRatio)
        
        # Define Seed and Ranges
        np.random.seed(r_seed)
        
        if NOISE:
            for i in range(len(realConcCC)):
                for j in range(len(realConcCC[i])):
                    realConcCC[i][j] = (realConcCC[i][j] + 
                              np.random.normal(0,ABS_NOISE_STD) +
                              np.random.normal(0,np.abs(realConcCC[i][j]*REL_NOISE_STD)))
            
            print("Control Coefficients with Noise Added")
            print(realConcCC)
        
        # Initualize Lists
        best_dist = []
        avg_dist = []
        med_dist = []
        top5_dist = []
        
        # Start Timing
        t1 = time.time()
        
        # Initialize
        ens_dist, ens_model, ens_rl, rl_track = initialize()
        
        dist_top_ind = np.argsort(ens_dist)
        dist_top = ens_dist[dist_top_ind]
        model_top = ens_model[dist_top_ind]
#        rl_top = ens_rl[dist_top_ind]
        
        print("Minimum distance: " + str(dist_top[0]))
        print("Average distance: " + str(np.average(dist_top)))
        best_dist.append(dist_top[0])
        avg_dist.append(np.average(dist_top))
        med_dist.append(np.median(dist_top))
        top5_dist.append(np.average(np.unique(dist_top)[:int(0.05*ens_size)]))
        
        breakFlag = False
        
        # TODO: Remove for loop
        for n in n_range:
            minind = np.argsort(ens_dist)[:pass_size]
            tarind = np.delete(np.arange(ens_size), minind)
            mut_p = 1/ens_dist[tarind]/np.sum(1/ens_dist[tarind])
            mut_ind = np.random.choice(tarind, size=mut_size-pass_size, 
                                               replace=False, p=mut_p)
            mut_ind = np.append(mut_ind, minind)
            mut_ind_inv = np.setdiff1d(np.arange(ens_size), mut_ind)
            
            eval_dist, eval_model, eval_rl = mutate_and_evaluate(ens_model[mut_ind], 
                                                                 ens_dist[mut_ind], 
                                                                 ens_rl[mut_ind])
            ens_model[mut_ind] = eval_model
            ens_dist[mut_ind] = eval_dist
            ens_rl[mut_ind] = eval_rl
            
    #        for tt in range(len(mut_ind)):
    #            r = te.loada(ens_model[mut_ind[tt]])
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
            
            rnd_dist, rnd_model, rnd_rl = random_gen(ens_model[mut_ind_inv], 
                                                     ens_dist[mut_ind_inv], 
                                                     ens_rl[mut_ind_inv])
            ens_model[mut_ind_inv] = rnd_model
            ens_dist[mut_ind_inv] = rnd_dist
            ens_rl[mut_ind_inv] = rnd_rl
            
            dist_top_ind = np.argsort(ens_dist)
            dist_top = ens_dist[dist_top_ind]
            model_top = ens_model[dist_top_ind]
#            rl_top = ens_rl[dist_top_ind]
            
            best_dist.append(dist_top[0])
            avg_dist.append(np.average(dist_top))
            med_dist.append(np.median(dist_top))
            top5_dist.append(np.average(np.unique(dist_top)[:int(0.05*ens_size)]))
            
            print("In generation: " + str(n + 1))
            print("Minimum distance: " + str(best_dist[-1]))
            print("Top 5 distance: " + str(top5_dist[-1]))
            print("Average distance: " + str(avg_dist[-1]))
            
            
    #        for tt in range(len(mut_ind_inv)):
    #            r = te.loada(ens_model[mut_ind_inv[tt]])
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
        if len(minInd[0]) == 0:
            minInd = np.array([[len(model_top) - 1]])
                
        if EXPORT_ALL_MODELS:
            minInd = np.array([[len(model_top) - 1]])
        
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
                pt.plotAllProgress([best_dist, avg_dist, med_dist, top5_dist], 
                                   labels=['Best', 'Avg', 'Median', 'Top 5 percent'],
                                   SAVE_PATH=os.path.join(EXPORT_PATH, 'images/AllConvergences.pdf'))
            else:
                pt.plotAllProgress([best_dist, avg_dist, med_dist, top5_dist], 
                                   labels=['Best', 'Avg', 'Median', 'Top 5 percent'])
            # TODO: Add polishing with fast optimizer 
            
            # Average residual
            if SAVE_PLOT:
                pt.plotResidual(realModel, ens_model, ens_dist, 
                                SAVE_PATH=os.path.join(EXPORT_PATH, 'images/average_residual.pdf'))
            else:
                pt.plotResidual(realModel, ens_model, ens_dist)
                
            # Distance histogram with KDE
            if SAVE_PLOT:
                pt.plotDistanceHistogramWithKDE(dist_top, log_dens, minInd, 
                                                SAVE_PATH=os.path.join(EXPORT_PATH, 
                                                                       'images/distance_hist_w_KDE.pdf'))
            else:
                pt.plotDistanceHistogramWithKDE(dist_top, log_dens, minInd)
                
            # RMSE histogram
#            r_real = te.loada(realModel)
#            k_real = r_real.getGlobalParameterValues()
#            
#            top_result_k = []
#            top_diff_k = []
#            
#            for i in ens_range:
#                r = te.loada(ens_model[np.argsort(ens_dist)[i]])
#                top_k = r.getGlobalParameterValues()
#                top_result_k.append(top_k)
#                try:
#                    top_diff_k.append(np.sqrt(np.divide(np.sum(np.square(np.subtract(
#                            k_real, top_k))),len(k_real))))
#                except:
#                    top_diff_k.append(np.sqrt(np.divide(np.sum(np.square(np.subtract(
#                            k_real, top_k[1:]))),len(k_real))))
#            
#            krmse = top_diff_k[:pass_size]
#            
#            plt.hist(krmse, bins=15, density=True)
#            plt.xlabel("RMSE", fontsize=15)
#            plt.ylabel("Normalized Frequency", fontsize=15)
#            plt.xticks(fontsize=15)
#            plt.yticks(fontsize=15)
#            if SAVE_PLOT:
#                plt.savefig(os.path.join(EXPORT_PATH, 'images/parameter_rmse_.pdf'), bbox_inches='tight')
#            plt.show()
            
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
                                      settings, t2-t1, rl_track, path=EXPORT_PATH)

        

