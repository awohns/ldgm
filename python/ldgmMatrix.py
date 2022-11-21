#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 10:18:06 2022

@author: psalehin
"""
import pandas as pd 
from scipy.sparse import csc_matrix
import numpy as np
from sksparse.cholmod import cholesky
from utility import *

class ldgmMatrix:
    

    def __init__(self, edge_list_path, snp_list_path):
        
        #load edge_list and snp_list - make sure they are for the same set of SNPs
        edge_list = pd.read_csv(edge_list_path, header =  None)
        snp_list = pd.read_csv(edge_list_path)
        
        assert(all(edge_list[0] <= edge_list[1]))
        # make matrix symmetric by adding (j,i,entry) for every edge (i,j,entry) where i<j
        edge_list2 = edge_list[[1,0,2]]
        edge_list2 = edge_list2[edge_list[0] < edge_list[1]]
        edge_list2.rename(columns = {1:0, 0:1}, inplace = True)
        edge_list = pd.concat((edge_list,edge_list2),axis=0)
        
        # compressed sparse column matrix
        self.matrix = csc_matrix( (edge_list[2],(edge_list[0], edge_list[1])) )
        
        # SNP list
        self.snps = snp_list
    
        self.name = edge_list_path
        
    def multiply(self, y, whichIndices):
        if np.ndim(y)==1:
            y = np.reshape(y,(y.shape[0],1))
         
        # handle indices not in precision matrix (i.e. all-zeros columns)
        Anz = self.matrix!=0
        v = Anz * np.ones(self.matrix.shape[0]) != 0
        otherIndices = np.logical_and (np.logical_not(whichIndices), v)
        
        # submatrices of A == self.matrix
        A_11 = sliceMatrix(self.matrix, otherIndices, otherIndices)
        factor = cholesky(A_11)
        A_00 = sliceMatrix(self.matrix, whichIndices, whichIndices) 
        A_10 = sliceMatrix(self.matrix, otherIndices, whichIndices)
        
        # x == (A/A_11) * y
        z = factor(A_10 * y)
        x = A_00 * y  - np.transpose(A_10) * z
        return x
  
    def divide(self, y, whichIndices):
        if np.ndim(y)==1:
            y = np.reshape(y,(y.shape[0],1))
         
        # handle indices not in precision matrix (i.e. all-zeros columns)
        Anz = self.matrix!=0
        incl = Anz * np.ones(self.matrix.shape[0]) != 0

        #diagonal elements should not be zero
        assert np.all(incl[whichIndices])
        
        # yp is y augmented with zeros           
        yp = np.zeros((self.matrix.shape[0],y.shape[1]), dtype=float, order='C')
        yp[whichIndices, :] = y

        #xp is x augmented with entries that can be ignored
        xp = np.zeros_like(yp)
        A_11 = sliceMatrix(self.matrix, incl, incl)
        factor = cholesky(A_11)

        xp[incl, :] = factor(yp[incl, :])
        x = xp[whichIndices, :]
        return x    
    