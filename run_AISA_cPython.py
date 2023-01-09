#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""----------AISA OPTIMIZATION-PYTHON WITH C-EXTENSION VERSION 1.0-------------

Created on Mon Aug 22 16:35:54 2022
@author: python-version:selamibeyhan
AISA opt. developer: Esref Bogar, Pamukkale University, Biomedical Eng.
Supervisor: Selami Beyhan, Izmir Demokrasi Univ. Electrical and Electronics Eng.
Cite method with paper: Bogar, Esref, and Selami Beyhan. "Adolescent Identity Search Algorithm
(AISA): A novel metaheuristic approach for solving optimization problems." Applied Soft 
Computing 95 (2020): 106503.

-------------------------------------------------------------------------------
To implement the AISA algorithm modify:
    i) cost function and stopping criteria
    ii) dimension and number of candidates
    iii) upper and lower bounds of parameters

Note: Developed code is much faster than pure python code when nPop>100.
Be careful about used memory size due to number of parameters and population size.
Note that there is no guarantee to find the global minimum!!
-------------------------------------------------------------------------------
"""
## Initial Parameters
import numpy as np
from numpy.random import rand
#import matplotlib.pyplot as plt
import time
import ChebyshevFeatureSelection
from Cost_Function import Cost_Function

#--------------------------------------MAIN------------------------------------
def main():
    #-----------------------------OPTIMIZATION LOOP----------------------------
    start_time = time.time()
    ub =  MaxValue*np.ones(Dim)     # ub: upper bound of parameters
    lb = -MaxValue*np.ones(Dim)     # lb: lower bound of parameters
    Candidates = np.zeros((nPop,Dim),dtype='d')
    Candidates_LSE = np.empty(shape=[0, Dim],dtype='d')
    X = np.zeros((nPop,Dim),dtype='d')
    canCost  = np.zeros(nPop,dtype='d')
    xcanCost = np.zeros(maxiter)
    ind_plt  = np.zeros(maxiter)
    cost_plt = np.zeros(maxiter)
    for i in range(0,nPop):
        Candidates[i] = np.multiply(rand(Dim),(ub-lb))+lb
        canCost[i] = Cost_Function(Candidates[i],Dim)
    flag = 0
    modify_term = np.zeros((1,Dim))
    for i in range(1,maxiter+1):
        if flag==0:
            #Following line calls the embedded C code
            Candidates_LSE = ChebyshevFeatureSelection.aisa_feature(nPop, Dim, Candidates.ravel().tolist(), canCost.ravel().tolist())
            ind1 = np.argsort(canCost)
            index = ind1[0]
        flag = 1
        for j in range(0,nPop):
            r1 = rand(1)
            if r1<=1/3:
                X[j] = Candidates[j]-(rand(1)*(Candidates[j]-Candidates_LSE))
            elif r1<=2/3 and r1>1/3:
                gg = np.random.randint(nPop,size=(1,1))
                if gg==index:
                    gg = np.random.randint(nPop,size=(1,1))
                X[j]= Candidates[j]-(rand(1))*(Candidates[gg]-Candidates[index])
            else:
                aa = np.random.randint(nPop,size=(1,1))
                bb = np.random.randint(Dim,size=(1,1))
                cc = np.random.randint(Dim,size=(1,1))
                modify_term = [Candidates[aa,bb]*np.ones(Dim)]-Candidates[j]
                X[j] = Candidates[j]+np.multiply(rand(1,Dim),modify_term)
            #Limit parameters--------------------------------------------------
            for k in range(0,Dim):
                if (X[j,k]>ub[k]) or (X[j,k]<lb[k]):
                    X[j,k] = lb[k] + rand(1)*(ub[k]-lb[k])
            xcanCost = Cost_Function(X[j],Dim)
            if xcanCost<canCost[j]:
                flag = 0
                Candidates[j] = X[j]
                canCost[j] = xcanCost
        ind2 = canCost.argmin()
        val2 = canCost[ind2]
        #--------Display the cost values at some batch-iterations----------------
        if np.remainder(i,maxiter/5)==0:
            print("Batch_number:",i,"Min_Cost:",val2)
        if val2<eps: break
        ind_plt[i-1] = i      # for plot indices
        cost_plt[i-1] = val2  # for plot values
    stop_time = time.time()
    # -------------------END OF OPTIMIZATION LOOP-------------------------------
    print('-----------------------------------------------')
    print('------------OPTIMIZATION RESULTS---------------')
    print('-----------------------------------------------')
    best_par = Candidates[ind2]
    print("Optimal parameters",best_par)
    print("Number of batch iterations:",i)
    print("Stopping criteria:",eps)
    print("Minimum cost value",val2)
    print("Elapsed time:",stop_time-start_time)

    #---------------------------------------------------------------------------
    #---------------------COST FUNCTION PLOT------------------------------------
    #plt.xlabel('Batch number')
    #plt.ylabel('Cost values')
    #plt.title('AISA optimization')
    #plt.plot(ind_plt[1:i-1], cost_plt[1:i-1],'-*b')
    #plt.show()

#-------------------------------------------------------------------------------
#------------------------MAIN FUNCTION------------------------------------------
if __name__ == "__main__":
    nPop    = 30                 # Number of populations
    maxiter = 100                # Number of batch iterations
    Dim = 2                      # Dimension of parameters
    MaxValue = 5                 # Max and min of values of parameters
    eps = 1e-8                   # Stopping criteria for cost value if needed.
    main()