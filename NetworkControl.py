# Import all necessary libraries
import pandas as pd
import numpy as np
import scipy
import scipy.stats as stats
from pathlib import Path
from scipy.io import loadmat
import scipy.sparse.linalg as sla
import scipy.linalg as la


def controllability_gramian(A, B, T = np.inf):
    '''
    Aim: Compute the causal controllability Gramian of the continuous time system.
    Input: 
        A the adjacency matrix as np.array
        B the control matrix (same size as the adjacency matrix specifying the control nodes as 1 and the non-control nodes as 0)
        T the horizon over which to compute the Gramian
    Output: 
        Matrix Wc: the inverse controllability gramian
    Details: 
    To compute the controllability Gramian for the whole matrix specify B as the identity matrix of A (as described in Cornblath et al. Nature Comms Biology 2020)
    Adjusted from Mark Muller's repository: https://github.com/markwmuller/controlpy
    '''
    
    # We need to solve the finite time Gramian
    # Boils down to solving an ODE:
    A = np.array(A,dtype=float)
    B = np.array(B,dtype=float)
    T = np.float(T)
    
    def gramian_ode(y, t0, A, B):
        temp = np.dot(scipy.linalg.expm(A*t0),B)
        dQ = np.dot(temp,np.conj(temp.T))
         
        return dQ.reshape((A.shape[0]**2,1))[:,0]
     
    y0 = np.zeros([A.shape[0]**2,1])[:,0]
    out = scipy.integrate.odeint(gramian_ode, y0, [0,T], args=(A,B))
    Q = out[1,:].reshape([A.shape[0], A.shape[0]])
    return Q


def minimalControlEnergy(adjacencyMatrix, WcI, x0, xf, T, normalise=True):
    ''' 
    Aim: Calculate minimal control energy to transition between states, using the whole brain as control matrix
    Inputs: 
        adjacencyMatrix: the structural matrix (N,N) as np.array
        WcI: the gramian inverse for control horizon T as a n (N,N) matrix
        x0: the initial state as a vector (N,1)
        xf: the final state as a vector (N,1) 
        T: the control time horizon    
    Output: 
        the minimal control energy needed to transition 
    Details: 
    Adapted from matlab code by Cornblath et al. Nature Comms Biology 2020, repository: https://github.com/ejcorn/brain_states/
    ''' 

    # Normalise if not already done
    if normalise==True:
        adjacencyMatrix = adjacencyMatrix / (np.full((adjacencyMatrix.shape), fill_value= sla.eigs(adjacencyMatrix)[0][0] + 1))
    # System size
    n = adjacencyMatrix.shape[0]
    # State transition to achieve
    exponentialMatrix = sla.expm(adjacencyMatrix*T)
    Phi = np.dot(exponentialMatrix, x0) - xf
    # Calculate energy needed to transition 
    E = np.sum(np.multiply(np.dot(WcI,Phi),(Phi)))
    return(E)

def normaliseMatrix(adjacencyMatrix):
    '''
    Aim: Normalise a connectivity matrix by maximum eigenvalue to make marginally stable
    Inputs:
        adjacencyMatrix: as an np.array (N*N nodes)
    Outputs:
        normalised matrix (N*N)
    '''
    normalisedMatrix = adjacencyMatrix / (np.full((adjacencyMatrix.shape), fill_value= sla.eigs(adjacencyMatrix)[0][0] + 1))
    return(normalisedMatrix)

def numberTransitions(affiliation, integratedLabel, segregatedLabel):
    '''
    Aim: calculate the number of transitions between integrated and segregated states 
        Note this is hard coded for the 2 states used in this approach
    Input: 
        affiliation: the list of labels for each of the fMRI volumes/windows as np.array
        integratedLabel: which label is the integrated state for the specific subject
        segregatedLabel: which label is the segregated state
    Output: 
        tuple (number of transitions Integrated to segregated, number of transitions Segregated to integrated)
    '''

    numTransitionsIntegratedtoSegregated = 0
    numTransitionsSegregatedtoIntegrated = 0

    for i in range(len(affiliation)-1):
        if (affiliation[i] == integratedLabel) & (affiliation[i+1] == segregatedLabel):
            numTransitionsIntegratedtoSegregated += 1 
        elif (affiliation[i] == segregatedLabel) & (affiliation[i+1] == integratedLabel):
            numTransitionsSegregatedtoIntegrated +=1
    return(numTransitionsIntegratedtoSegregated, numTransitionsSegregatedtoIntegrated)


def dwellTime(affiliation, integratedLabel, segregatedLabel):   
    '''
    Aim: calculate the dwell time in each of the 2 states (integrated and segregated)
    Input: 
        affiliation: the list of labels for each of the fMRI volumes/windows as np.array
        integratedLabel: which label is the integrated state for the specific subject
        segregatedLabel: which label is the segregated state
    Output: 
        tuple (mean dwell time in integrated, mean dwell time in segregated state)
    '''
    prevNumber =  affiliation[0]
    dwellTimeIntegrated = []
    dwellTimeSegregated = []
    dwellInt = 1
    dwellSeg = 1

    for i in range(1, len(affiliation)):
        if affiliation[i] == prevNumber:
            if prevNumber == integratedLabel: 
                dwellInt += 1
            elif prevNumber == segregatedLabel:
                dwellSeg +=1
        elif affiliation[i] != prevNumber:
            if prevNumber == integratedLabel:
                dwellTimeIntegrated.append(dwellInt)
                dwellInt = 1
                prevNumber = segregatedLabel
            elif prevNumber == segregatedLabel:
                dwellTimeSegregated.append(dwellSeg)
                dwellSeg = 1
                prevNumber = integratedLabel
    dwellIntegrated = np.mean(dwellTimeIntegrated)
    dwellSegregated = np.mean(dwellTimeSegregated)
    return(dwellIntegrated, dwellSegregated) 