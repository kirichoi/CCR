# -*- coding: utf-8 -*-

# Generate random mass-action networks using uniuni, biuni, unibi and bibi
# The following reaction patterns are currently not allowed:
# X -> X

# X + Y -> X
# Y + X -> X 
# X + X -> X

# X -> X + Y
# X -> Y + X
# X -> X + X

# X + X -> X + X
# X + Y  -> X + Y
# Y + X -> X + Y
# Y + X -> Y + X
# X + Y -> X + Y
# X + Y -> Y + X

# How to use:
#   rl =  generateReactionList (6, 5)
#   st = getFullStoichiometryMatrix (rl)
#   stt = removeBoundaryNodes (st)
#   if len (stt[1]) > 0:
#      antStr = genAntimonyScript (stt[1], stt[2], rl)

# floating and boundary Ids are represented as integers

import tellurium as te
import random
import numpy as np
import copy
import analysis

class TReactionType:
    UNIUNI_N = 0
    UNIUNI_A = 1
    UNIUNI_I = 2
    UNIUNI_AI = 3
    BIUNI_N = 4
    BIUNI_A = 5
    BIUNI_I = 6
    BIUNI_AI = 7
    UNIBI_N = 8
    UNIBI_A = 9
    UNIBI_I = 10
    UNIBI_AI = 11
    BIBI_N = 12
    BIBI_A = 13
    BIBI_I = 14
    BIBI_AI = 15
      
    
class TReactionProbabilities:
    UniUni = 0.7
    BiUni = 0.125
    UniBi = 0.125
    BiBI  = 0.05
    
    
class RateLawProb:
    default = 0.83
    inhib = 0.08
    activ = 0.08
    inhibactiv = 0.01

#def pickRateLawType():
#    rt = np.random.random()
#    if rt < RateLawProb.default:
#        return 0
#    elif rt < RateLawProb.default + RateLawProb.inhib:
#        return 1
#    elif rt < RateLawProb.default + RateLawProb.inhib + RateLawProb.activ:
#        return 2
#    return 3
     

def pickReactionType():
    rt1 = np.random.random()
    rt2 = np.random.random()
    if rt1 < TReactionProbabilities.UniUni:
        if rt2 < RateLawProb.default:
            return TReactionType.UNIUNI_N
        elif rt2 < RateLawProb.default + RateLawProb.activ:
            return TReactionType.UNIUNI_A
        elif rt2 < RateLawProb.default + RateLawProb.activ + RateLawProb.inhib:
            return TReactionType.UNIUNI_I
        else:
            return TReactionType.UNIUNI_AI
    elif rt1 < TReactionProbabilities.UniUni + TReactionProbabilities.BiUni:
        if rt2 < RateLawProb.default:
            return TReactionType.BIUNI_N
        elif rt2 < RateLawProb.default + RateLawProb.activ:
            return TReactionType.BIUNI_A
        elif rt2 < RateLawProb.default + RateLawProb.activ + RateLawProb.inhib:
            return TReactionType.BIUNI_I
        else:
            return TReactionType.BIUNI_AI
    elif rt1 < TReactionProbabilities.UniUni + TReactionProbabilities.BiUni + TReactionProbabilities.UniBi:
        if rt2 < RateLawProb.default:
            return TReactionType.UNIBI_N
        elif rt2 < RateLawProb.default + RateLawProb.activ:
            return TReactionType.UNIBI_A
        elif rt2 < RateLawProb.default + RateLawProb.activ + RateLawProb.inhib:
            return TReactionType.UNIBI_I
        else:
            return TReactionType.UNIBI_AI
    else:
        if rt2 < RateLawProb.default:
            return TReactionType.BIBI_N
        elif rt2 < RateLawProb.default + RateLawProb.activ:
            return TReactionType.BIBI_A
        elif rt2 < RateLawProb.default + RateLawProb.activ + RateLawProb.inhib:
            return TReactionType.BIBI_I
        else:
            return TReactionType.BIBI_AI


# Generates a reaction network in the form of a reaction list
# reactionList = [nSpecies, reaction, ....]
# reaction = [reactionType, [list of reactants], [list of product], rateConstant]
def generateReactionList(nSpecies, nReactions, boundaryIdx):
    connected = False
    
    while not connected:
        reactionList = []
        for r in range(nReactions):
            rct = [col[1] for col in reactionList]
            prd = [col[2] for col in reactionList]
            #rateConstant = 0.5  #np.random.uniform(1e-3, 1.)
            rt = pickReactionType()
            
            if rt <= TReactionType.UNIUNI_AI:
                # UniUni
                rct_id = np.random.choice(np.arange(nSpecies), size=1)
                prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=1)
                all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                
                while ((np.any(np.isin(rct_id, boundaryIdx))) and (np.any(np.isin(prd_id, boundaryIdx)))) or (len(set(all_rct) & set(all_prd)) > 0):
                    rct_id = np.random.choice(np.arange(nSpecies), size=1)
                    prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=1)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                    all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                    
                    while len(set(all_rct) & set(all_prd)) > 0:
                        rct_id = np.random.choice(np.arange(nSpecies), size=1)
                        prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=1)
                        all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                        all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                
                if rt == TReactionType.UNIUNI_N:
                    act_id = []
                    inhib_id = []
                elif rt == TReactionType.UNIUNI_A:
                    act_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=1).tolist()
                    inhib_id = []
                    if len(act_id) == 0:
                        rt = TReactionType.UNIUNI_N
                elif rt == TReactionType.UNIUNI_I:
                    act_id = []
                    inhib_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=1).tolist()
                    if len(inhib_id) == 0:
                        rt = TReactionType.UNIUNI_N
                else:
                    reg_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=2)
                    act_id = [reg_id[0]]
                    inhib_id = [reg_id[1]]
                    if len(reg_id) == 0:
                        rt = TReactionType.UNIUNI_N
                
                reactionList.append ([rt, [rct_id[0]], [prd_id[0]], act_id, inhib_id])
    
            elif (rt > TReactionType.UNIUNI_AI) and (rt <= TReactionType.BIUNI_AI):
                # BiUni
                rct_id = np.random.choice(np.arange(nSpecies), size=2, replace=True)
                prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=1)
                all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                
                while ((np.any(np.isin(rct_id, boundaryIdx))) and (np.any(np.isin(prd_id, boundaryIdx)))) or (len(set(all_rct) & set(all_prd)) > 0):
                    rct_id = np.random.choice(np.arange(nSpecies), size=2, replace=True)
                    # pick a product but only products that don't include the reactants
                    prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=1)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                    all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                    
                    while len(set(all_rct) & set(all_prd)) > 0:
                        rct_id = np.random.choice(np.arange(nSpecies), size=2, replace=True)
                        prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=1)
                        all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                        all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0]]]
                
                if rt == TReactionType.BIUNI_N:
                    act_id = []
                    inhib_id = []
                elif rt == TReactionType.BIUNI_A:
                    act_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=1).tolist()
                    inhib_id = []
                    if len(act_id) == 0:
                        rt = TReactionType.BIUNI_N
                elif rt == TReactionType.BIUNI_I:
                    act_id = []
                    inhib_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=1).tolist()
                    if len(inhib_id) == 0:
                        rt = TReactionType.BIUNI_N
                else:
                    reg_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=2)
                    act_id = [reg_id[0]]
                    inhib_id = [reg_id[1]]
                    if len(reg_id) == 0:
                        rt = TReactionType.BIUNI_N
                
                reactionList.append ([rt, [rct_id[0], rct_id[1]], [prd_id[0]], act_id, inhib_id]) 
            
            elif (rt > TReactionType.BIUNI_AI) and (rt <= TReactionType.UNIBI_AI):
                # UniBi
                rct_id = np.random.choice(np.arange(nSpecies), size=1)
                prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=2, replace=True)
                all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                
                while ((np.any(np.isin(rct_id, boundaryIdx))) and (np.any(np.isin(prd_id, boundaryIdx)))) or (len(set(all_rct) & set(all_prd)) > 0):
                    rct_id = np.random.choice(np.arange(nSpecies), size=1)
                    # pick a product but only products that don't include the reactant
                    prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=2, replace=True)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                    all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                    
                    while len(set(all_rct) & set(all_prd)) > 0:
                        rct_id = np.random.choice(np.arange(nSpecies), size=1)
                        prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=2, replace=True)
                        all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0]]]
                        all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                        
                if rt == TReactionType.UNIBI_N:
                    act_id = []
                    inhib_id = []
                elif rt == TReactionType.UNIBI_A:
                    act_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=1).tolist()
                    inhib_id = []
                    if len(act_id) == 0:
                        rt = TReactionType.UNIBI_N
                elif rt == TReactionType.UNIBI_I:
                    act_id = []
                    inhib_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=1).tolist()
                    if len(inhib_id) == 0:
                        rt = TReactionType.UNIBI_N
                else:
                    reg_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=2)
                    act_id = [reg_id[0]]
                    inhib_id = [reg_id[1]]
                    if len(reg_id) == 0:
                        rt = TReactionType.UNIBI_N
                        
                reactionList.append ([rt, [rct_id[0]], [prd_id[0], prd_id[1]], act_id, inhib_id]) 
            
            else:
                # BiBi
                rct_id = np.random.choice(np.arange(nSpecies), size=2, replace=True)
                prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=2, replace=True)
                all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                
                while ((np.any(np.isin(rct_id, boundaryIdx))) and (np.any(np.isin(prd_id, boundaryIdx)))) or (len(set(all_rct) & set(all_prd)) > 0):
                    rct_id = np.random.choice(np.arange(nSpecies), size=2, replace=True)
                    # pick a product but only products that don't include the reactant
                    prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=2, replace=True)
                    
                    # Search for potentially identical reactions
                    all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                    all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                    
                    while len(set(all_rct) & set(all_prd)) > 0:
                        rct_id = np.random.choice(np.arange(nSpecies), size=2, replace=True)
                        prd_id = np.random.choice(np.delete(np.arange(nSpecies), rct_id), size=2, replace=True)
                        all_rct = [i for i, x in enumerate(rct) if x == [rct_id[0],
                                                       rct_id[1]] or x == [rct_id[1], rct_id[0]]]
                        all_prd = [i for i, x in enumerate(prd) if x == [prd_id[0],
                                                       prd_id[1]] or x == [prd_id[1], prd_id[0]]]
                
                if rt == TReactionType.BIBI_N:
                    act_id = []
                    inhib_id = []
                elif rt == TReactionType.BIBI_A:
                    act_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=1).tolist()
                    inhib_id = []
                    if len(act_id) == 0:
                        rt = TReactionType.BIBI_N
                elif rt == TReactionType.BIBI_I:
                    act_id = []
                    inhib_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=1).tolist()
                    if len(inhib_id) == 0:
                        rt = TReactionType.BIBI_N
                else:
                    reg_id = np.random.choice(np.delete(np.delete(np.arange(nSpecies), rct_id), prd_id), size=2)
                    act_id = [reg_id[0]]
                    inhib_id = [reg_id[1]]
                    if len(reg_id) == 0:
                        rt = TReactionType.UNIBI_N
                
                reactionList.append ([rt, [rct_id[0], rct_id[1]], [prd_id[0], prd_id[1]], act_id, inhib_id])
        
        connected = analysis.isConnected(reactionList)
        
    return reactionList
    

# Include boundary and floating species
# Returns a list:
# [New Stoichiometry matrix, list of floatingIds, list of boundaryIds]
def getFullStoichiometryMatrix(reactionList, ns):
    reactionListCopy = copy.deepcopy(reactionList)
    st = np.zeros((ns, len(reactionListCopy)))
    
    for index, rind in enumerate(reactionListCopy):
        if (rind[0] <= TReactionType.UNIUNI_AI):
            # UniUni
            reactant = reactionListCopy[index][1][0]
            st[reactant, index] = st[reactant, index] - 1
            product = reactionListCopy[index][2][0]
            st[product, index] = st[product, index] + 1
     
        elif (rind[0] > TReactionType.UNIUNI_AI) and (rind[0] <= TReactionType.BIUNI_AI):
            # BiUni
            reactant1 = reactionListCopy[index][1][0]
            st[reactant1, index] = st[reactant1, index] - 1
            reactant2 = reactionListCopy[index][1][1]
            st[reactant2, index] = st[reactant2, index] - 1
            product = reactionListCopy[index][2][0]
            st[product, index] = st[product, index] + 1

        elif (rind[0] > TReactionType.BIUNI_AI) and (rind[0] <= TReactionType.UNIBI_AI):
            # UniBi
            reactant1 = reactionListCopy[index][1][0]
            st[reactant1, index] = st[reactant1, index] - 1
            product1 = reactionListCopy[index][2][0]
            st[product1, index] = st[product1, index] + 1
            product2 = reactionListCopy[index][2][1]
            st[product2, index] = st[product2, index] + 1

        else:
            # BiBi
            reactant1 = reactionListCopy[index][1][0]
            st[reactant1, index] = st[reactant1, index] - 1
            reactant2 = reactionListCopy[index][1][1]
            st[reactant2, index] = st[reactant2, index] - 1
            product1 = reactionListCopy[index][2][0]
            st[product1, index] = st[product1, index] + 1
            product2 = reactionListCopy[index][2][1]
            st[product2, index] = st[product2, index] + 1

    return st
        

# Removes boundary or orphan species from stoichiometry matrix
def removeBoundaryNodes(st):
    
    dims = st.shape
    
    nSpecies = dims[0]
    nReactions = dims[1]
    
    speciesIds = np.arange (nSpecies)
    indexes = []
    orphanSpecies = []
    countBoundarySpecies = 0
    for r in range(nSpecies): 
        # Scan across the columns, count + and - coefficients
        plusCoeff = 0; minusCoeff = 0
        for c in range(nReactions):
            if st[r,c] < 0:
                minusCoeff = minusCoeff + 1
            elif st[r,c] > 0:
                plusCoeff = plusCoeff + 1
        if plusCoeff == 0 and minusCoeff == 0:
            # No reaction attached to this species
            orphanSpecies.append (r)
        elif plusCoeff == 0 and minusCoeff != 0:
            # Species is a source
            indexes.append (r)
            countBoundarySpecies = countBoundarySpecies + 1
        elif minusCoeff == 0 and plusCoeff != 0:
            # Species is a sink
            indexes.append (r)
            countBoundarySpecies = countBoundarySpecies + 1

    floatingIds = np.delete(speciesIds, indexes+orphanSpecies, axis=0).astype('int')
    floatingIds = floatingIds.tolist()

    boundaryIds = indexes
    return [np.delete(st, indexes + orphanSpecies, axis=0), floatingIds, boundaryIds]


def generateRateLaw(rl, floatingIds, boundaryIds, rlt, Jind):
    
    Klist = []
    
    T = ''
    D = ''
    Rreg = ''
    Dreg = ''
    
    # T
    T = T + '(Kf_' + str(Jind) + '*'
    Klist.append('Kf_' + str(Jind))
    Klist.append('h_' + str(Jind))
    
    T = T + '('
    for i in range(len(rl[Jind][1])):
        T = T + '(S' + str(rl[Jind][1][i]) + '/Km_' + str(rl[Jind][1][i]) + ')^h_' + str(Jind)
        Klist.append('Km_' + str(rl[Jind][1][i]))
        if i < len(rl[Jind][1]) - 1:
            T = T + '*'
    T = T + ')'
    
    T = T + '-Kr_' + str(Jind) + '*'
    Klist.append('Kr_' + str(Jind))
    
    T = T + '('
    for i in range(len(rl[Jind][2])):
        T = T + '(S' + str(rl[Jind][2][i]) + '/Km_' + str(rl[Jind][2][i]) +')^h_' + str(Jind)
        Klist.append('Km_' + str(rl[Jind][2][i]))
        if i < len(rl[Jind][2]) - 1:
            T = T + '*'
            
    T = T + '))'
        
    # D
    D = D + '('
    
    for i in range(len(rl[Jind][1])):
        D = D + '((1 + (S' + str(rl[Jind][1][i]) + '/Km_' + str(rl[Jind][1][i]) + '))^h_' + str(Jind) + ')'
        Klist.append('Km_' + str(rl[Jind][1][i]))
        if i < len(rl[Jind][1]) - 1:
            D = D + '*'
    
    D = D + '+'
    
    for i in range(len(rl[Jind][2])):
        D = D + '((1 + (S' + str(rl[Jind][2][i]) + '/Km_' + str(rl[Jind][2][i]) + '))^h_' + str(Jind) + ')'
        Klist.append('Km_' + str(rl[Jind][2][i]))
        if i < len(rl[Jind][2]) - 1:
            D = D + '*'
    
    D = D + '-1)'
        
    #Rreg
    if (rlt == 1) or (rlt == 3):
        pass
    
    #Dreg
    if (rlt == 2) or (rlt == 3):
        pass
    
    
    rateLaw = Rreg + T + '/(' + D +  Dreg + ')'
        
    return rateLaw, Klist


def generateSimpleRateLaw(rl, floatingIds, boundaryIds, Jind):
    
    Klist = []
    
    T = ''
    D = ''
    ACT = ''
    INH = ''
    
    # T
    T = T + '(Kf_' + str(Jind) + '*'
    Klist.append('Kf_' + str(Jind))
    
    for i in range(len(rl[Jind][1])):
        T = T + 'S' + str(rl[Jind][1][i])
        if i < len(rl[Jind][1]) - 1:
            T = T + '*'
    
    T = T + ' - Kr_' + str(Jind) + '*'
    Klist.append('Kr_' + str(Jind))
    
    for i in range(len(rl[Jind][2])):
        T = T + 'S' + str(rl[Jind][2][i])
        if i < len(rl[Jind][2]) - 1:
            T = T + '*'
            
    T = T + ')'
        
    # D
    D = D + '1 + '
    
    for i in range(len(rl[Jind][1])):
        D = D + 'S' + str(rl[Jind][1][i])
        if i < len(rl[Jind][1]) - 1:
            D = D + '*'
    
    D = D + ' + '
    
    for i in range(len(rl[Jind][2])):
        D = D + 'S' + str(rl[Jind][2][i])
        if i < len(rl[Jind][2]) - 1:
            D = D + '*'
    
    # Activation
    if (len(rl[Jind][3]) > 0):
        ACT = ACT + ' + '
        for i in range(len(rl[Jind][3])):
            ACT = ACT + '1/S' + str(rl[Jind][3][i])
            if i < len(rl[Jind][3]) - 1:
                ACT = ACT + ' + '
    
    # Inhibition
    if (len(rl[Jind][4]) > 0):
        INH = INH + ' + '
        for i in range(len(rl[Jind][4])):
            INH = INH + 'S' + str(rl[Jind][4][i])
            if i < len(rl[Jind][4]) - 1:
                INH = INH + ' + '
    
    
    rateLaw = T + '/(' + D + ACT + INH + ')'
        
    return rateLaw, Klist


def generateAntimony(floatingIds, boundaryIds, stt1, stt2, reactionList, boundary_init=None):
    reactionListCopy = copy.deepcopy(reactionList)
    Klist = []
    
    real = np.append(floatingIds, boundaryIds)
    if type(real[0]) == 'str' or type(real[0]) == np.str_:
        real = [s.strip('S') for s in real]
    real = list(map(int, real))
    tar = stt1 + stt2
    tar = list(map(int, tar))
    
    # List species
    antStr = ''
    if len (floatingIds) > 0:
        antStr = antStr + 'var ' + str(floatingIds[0])
        for index in floatingIds[1:]:
            antStr = antStr + ', ' + str(index)
        antStr = antStr + ';\n'
    
    if len (boundaryIds) > 0:
        antStr = antStr + 'const ' + str(boundaryIds[0])
        for index in boundaryIds[1:]:
            antStr = antStr + ', ' + str(index)
        antStr = antStr + ';\n'

    # List reactions
    for index, rind in enumerate(reactionListCopy):
        if (rind[0] <= TReactionType.UNIUNI_AI):
            # UniUni
            antStr = antStr + 'J' + str(index) + ': S' + str(real[tar.index(reactionListCopy[index][1][0])])
            antStr = antStr + ' -> '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][0])])
            antStr = antStr + '; '#k' + str(index) + '*S' + str(real[tar.index(reactionListCopy[index][1][0])])
            RateLaw, klist_i = generateRateLaw(reactionList, floatingIds, boundaryIds, 0, index)
            antStr = antStr + RateLaw
            Klist.append(klist_i)
        elif (rind[0] > TReactionType.UNIUNI_AI) and (rind[0] <= TReactionType.BIUNI_AI):
            # BiUni
            antStr = antStr + 'J' + str(index) + ': S' + str(real[tar.index(reactionListCopy[index][1][0])])
            antStr = antStr + ' + '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][1][1])])
            antStr = antStr + ' -> '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][0])])
            antStr = antStr + '; '#k' + str(index) + '*S' + str(real[tar.index(reactionListCopy[index][1][0])]) + '*S' + str(real[tar.index(reactionListCopy[index][1][1])])
            RateLaw, klist_i = generateRateLaw(reactionList, floatingIds, boundaryIds, 0, index)
            antStr = antStr + RateLaw
            Klist.append(klist_i)
        elif (rind[0] > TReactionType.BIIUNI_AI) and (rind[0] <= TReactionType.UNIBI_AI):
            # UniBi
            antStr = antStr + 'J' + str(index) + ': S' + str(real[tar.index(reactionListCopy[index][1][0])])
            antStr = antStr + ' -> '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][0])])
            antStr = antStr + ' + '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][1])])
            antStr = antStr + '; '#k' + str(index) + '*S' + str(real[tar.index(reactionListCopy[index][1][0])])
            RateLaw, klist_i = generateRateLaw(reactionList, floatingIds, boundaryIds, 0, index)
            antStr = antStr + RateLaw
            Klist.append(klist_i)
        else:
            # BiBi
            antStr = antStr + 'J' + str(index) + ': S' + str(real[tar.index(reactionListCopy[index][1][0])])
            antStr = antStr + ' + '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][1][1])])
            antStr = antStr + ' -> '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][0])])
            antStr = antStr + ' + '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][1])])
            antStr = antStr + '; '#k' + str(index) + '*S' + str(real[tar.index(reactionListCopy[index][1][0])]) + '*S' + str(real[tar.index(reactionListCopy[index][1][1])])
            RateLaw, klist_i = generateRateLaw(reactionList, floatingIds, boundaryIds, 0, index)
            antStr = antStr + RateLaw
            Klist.append(klist_i)
        antStr = antStr + ';\n'

    # List rate constants
    antStr = antStr + '\n'
    Klist_f = [item for sublist in Klist for item in sublist]
    Klist_f = np.unique(Klist_f)
    for i in range(len(Klist_f)):
        if Klist_f[i].startswith('Km'):
            antStr = antStr + Klist_f[i] + ' = 1\n'
        elif Klist_f[i].startswith('h'):
            antStr = antStr + Klist_f[i] + ' = 1\n'
        elif Klist_f[i].startswith('Kf'):
            antStr = antStr + Klist_f[i] + ' = 1\n'
        elif Klist_f[i].startswith('Kr'):
            antStr = antStr + Klist_f[i] + ' = 0.5\n'
        
    # Initialize boundary species
    antStr = antStr + '\n'
    if type(boundary_init) == type(None):
        for index, bind in enumerate(boundaryIds):
            antStr = antStr + str(bind) + ' = ' + str(np.random.randint (1,6)) + '\n'
    else:
        for index, bind in enumerate(boundaryIds):
            antStr = antStr + str(bind) + ' = ' + str(boundary_init[index]) + '\n'
    
    # Initialize floating species
    for index, find in enumerate(floatingIds):
        antStr = antStr + str(find) + ' = ' + '0\n'
        
    return antStr
     

def generateParameterBoundary(glgp):
    
    pBound = []
    
    for i in range(len(glgp)):
        if glgp[i].startswith('Km'):
            pBound.append((1e-3, 10.))
        elif glgp[i].startswith('h'):
            pBound.append((1e-1, 4.))
        elif glgp[i].startswith('Kf'):
            pBound.append((1e-3, 10.))
        elif glgp[i].startswith('Kr'):
            pBound.append((1e-3, 10.))

    return pBound
    

def generateLinearChainAnt(ns):
    order = np.random.sample(range(ns), ns)
    
    antStr = ''
    antStr = antStr + 'var S' + str(order[1])
    
    for i in range(ns - 3):
        antStr = antStr + ', S' + str(order[i + 2])
    
    antStr = antStr + '\n'
    antStr = antStr + 'const S' + str(order[0]) + ', S' + str(order[-1])
    antStr = antStr + ';\n'
    
    for i in range(ns - 1):
        antStr = antStr + 'S' + str(order[i]) + ' -> S' + str(order[i + 1]) + '; k' + str(order[i]) + '*S' + str(order[i]) + '\n'
    
    antStr = antStr + '\n'
    antStr = antStr + 'S' + str(order[0]) + ' = ' + str(random.randint (1,6)) + ';\n'
    
    for i in range(ns - 1):
        antStr = antStr + 'k' + str(order[i]) + ' = ' + str(random.random()) + ';\n'
        
    return antStr

    
