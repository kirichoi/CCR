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

class TReactionType:
    UNIUNI = 0
    BIUNI = 1
    UNIBI = 2
    BIBI = 3
      
class TReactionProbabilities:
    UniUni = 0.7
    BiUni = 0.125
    UniBi = 0.125
    BiBI  = 0.05
     
      
def pickReactionType():
    rt = np.random.random()
    if rt < TReactionProbabilities.UniUni:
        return TReactionType.UNIUNI
    elif rt < TReactionProbabilities.UniUni + TReactionProbabilities.BiUni:
        return TReactionType.BIUNI
    elif rt < TReactionProbabilities.UniUni + TReactionProbabilities.BiUni + TReactionProbabilities.UniBi:
        return TReactionType.UNIBI
    return TReactionType.BIBI
    

# Generates a reaction network in the form of a reaction list
# reactionList = [nSpecies, reaction, ....]
# reaction = [reactionType, [list of reactants], [list of product], rateConstant]
def generateReactionList(nSpecies, nReactions):
    reactionList = []
    for r in range(nReactions):
        rct = [col[1] for col in reactionList]
        prd = [col[2] for col in reactionList]
        rateConstant = 0.5#np.random.uniform(1e-3, 1.)
        rt = pickReactionType()
        if rt == TReactionType.UNIUNI:
            # UniUni
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
            
            reactionList.append ([rt, [rct_id[0]], [prd_id[0]], rateConstant]) 

        elif rt == TReactionType.BIUNI:
            # BiUni
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
            
            reactionList.append ([rt, [rct_id[0], rct_id[1]], [prd_id[0]], rateConstant]) 
        
        elif rt == TReactionType.UNIBI:
            # UniBi
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
            
            reactionList.append ([rt, [rct_id[0]], [prd_id[0], prd_id[1]], rateConstant]) 
        
        else:
            # BiBi
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
            
            reactionList.append ([rt, [rct_id[0], rct_id[1]], [prd_id[0], prd_id[1]], rateConstant])

    return reactionList
    

# Include boundary and floating species
# Returns a list:
# [New Stoichiometry matrix, list of floatingIds, list of boundaryIds]
def getFullStoichiometryMatrix(reactionList, ns):
    reactionListCopy = copy.deepcopy(reactionList)
    st = np.zeros((ns, len(reactionListCopy)))
    
    for index, rind in enumerate(reactionListCopy):
        if rind[0] == TReactionType.UNIUNI:
            # UniUni
            reactant = reactionListCopy[index][1][0]
            st[reactant, index] = st[reactant, index] - 1
            product = reactionListCopy[index][2][0]
            st[product, index] = st[product, index] + 1
     
        elif rind[0] == TReactionType.BIUNI:
            # BiUni
            reactant1 = reactionListCopy[index][1][0]
            st[reactant1, index] = st[reactant1, index] - 1
            reactant2 = reactionListCopy[index][1][1]
            st[reactant2, index] = st[reactant2, index] - 1
            product = reactionListCopy[index][2][0]
            st[product, index] = st[product, index] + 1

        elif rind[0] == TReactionType.UNIBI:
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


def generateRateLaw(rt, rl):
    pass


def generateAntimony(floatingIds, boundaryIds, stt1, stt2, reactionList, boundary_init=None):
    reactionListCopy = copy.deepcopy(reactionList)
    
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
        if rind[0] == TReactionType.UNIUNI:
            # UniUni
            antStr = antStr + 'J' + str(index) + ': S' + str(real[tar.index(reactionListCopy[index][1][0])])
            antStr = antStr + ' -> '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][0])])
            antStr = antStr + '; k' + str(index) + '*S' + str(real[tar.index(reactionListCopy[index][1][0])])
        elif rind[0] == TReactionType.BIUNI:
            # BiUni
            antStr = antStr + 'J' + str(index) + ': S' + str(real[tar.index(reactionListCopy[index][1][0])])
            antStr = antStr + ' + '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][1][1])])
            antStr = antStr + ' -> '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][0])])
            antStr = antStr + '; k' + str(index) + '*S' + str(real[tar.index(reactionListCopy[index][1][0])]) + '*S' + str(real[tar.index(reactionListCopy[index][1][1])])
        elif rind[0] == TReactionType.UNIBI:
            # UniBi
            antStr = antStr + 'J' + str(index) + ': S' + str(real[tar.index(reactionListCopy[index][1][0])])
            antStr = antStr + ' -> '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][0])])
            antStr = antStr + ' + '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][1])])
            antStr = antStr + '; k' + str(index) + '*S' + str(real[tar.index(reactionListCopy[index][1][0])])
        else:
            # BiBi
            antStr = antStr + 'J' + str(index) + ': S' + str(real[tar.index(reactionListCopy[index][1][0])])
            antStr = antStr + ' + '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][1][1])])
            antStr = antStr + ' -> '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][0])])
            antStr = antStr + ' + '
            antStr = antStr + 'S' + str(real[tar.index(reactionListCopy[index][2][1])])
            antStr = antStr + '; k' + str(index) + '*S' + str(real[tar.index(reactionListCopy[index][1][0])]) + '*S' + str(real[tar.index(reactionListCopy[index][1][1])])
        antStr = antStr + ';\n'

    # List rate constants
    antStr = antStr + '\n'
    for index, rind in enumerate(reactionListCopy):
        antStr = antStr + 'k' + str(index) + ' = ' + str(rind[3]) + '\n'
        
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

    
