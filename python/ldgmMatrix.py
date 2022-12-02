#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import pandas as pd 
from scipy.sparse import csc_matrix
import numpy as np
from sksparse.cholmod import cholesky
from utility import *

class ldgmMatrix:
    

    def __init__(self, edge_list_path = None, snp_list_path = None):
        
        self.matrix = None
        self.name = None
        self.snps = None
        self.nz_index = None
        
        if edge_list_path is not None:
            
            #load edge_list and snp_list - make sure they are for the same set of SNPs
            edge_list = pd.read_csv(edge_list_path, header =  None)
            
            assert(all(edge_list[0] <= edge_list[1]))
            # make matrix symmetric by adding (j,i,entry) for every edge (i,j,entry) where i<j
            edge_list2 = edge_list[[1,0,2]]
            edge_list2 = edge_list2[edge_list[0] < edge_list[1]]
            edge_list2.rename(columns = {1:0, 0:1}, inplace = True)
            edge_list = pd.concat((edge_list,edge_list2),axis=0)
            
            # compressed sparse column matrix
            self.matrix = csc_matrix( (edge_list[2],(edge_list[0], edge_list[1])) )
                    
            self.name = edge_list_path
            
            Anz = self.matrix!=0
            self.nz_index = Anz * np.ones(self.matrix.shape[0]) != 0
        
        if snp_list_path is not None:
            # SNP list
            snp_list = pd.read_csv(snp_list_path)
            assert "index" in snp_list.columns, "SNP list should have a column named 'index'"
            self.snps = snp_list
        
    def multiply(self, y, whichIndices):
        if np.ndim(y)==1:
            y = np.reshape(y,(y.shape[0],1))
        
        assert all(self.nz_index(whichIndices)), "Matrix should have nonzero diagonal entries for all nonmissing SNPs"
        # handle indices not in precision matrix (i.e. all-zeros columns)
        otherIndices = np.logical_and (np.logical_not(whichIndices), self.nz_idnex)
        
        # submatrices of A == self.matrix
        A_00 = sliceMatrix(self.matrix, otherIndices, otherIndices)
        factor = cholesky(A_00)
        A_11 = sliceMatrix(self.matrix, whichIndices, whichIndices) 
        A_01 = sliceMatrix(self.matrix, otherIndices, whichIndices)
        
        # x == (A/A_00) * y
        z = factor(A_01 * y)
        x = A_11 * y  - np.transpose(A_01) * z
        return x
  
    def divide(self, y, whichIndices):
        assert all(self.nz_index(whichIndices)), "Matrix should have nonzero diagonal entries for all nonmissing SNPs"
        
        if np.ndim(y)==1:
            y = np.reshape(y,(y.shape[0],1))
         
        #diagonal elements should not be zero
        assert np.all(self.nz_idnex[whichIndices])
        
        # yp is y augmented with zeros           
        yp = np.zeros((self.matrix.shape[0],y.shape[1]), dtype=float, order='C')
        yp[whichIndices, :] = y

        #xp is x augmented with entries that can be ignored
        xp = np.zeros_like(yp)
        A_11 = sliceMatrix(self.matrix, self.nz_idnex, self.nz_idnex)
        factor = cholesky(A_11)

        xp[self.nz_idnex, :] = factor(yp[self.nz_idnex, :])
        x = xp[whichIndices, :]
        return x

    def root_divide(self, y, whichIndices):
        assert all(self.nz_index[whichIndices]), "Matrix should have nonzero diagonal entries for all nonmissing SNPs"
        
        if np.ndim(y)==1:
            y = np.reshape(y,(y.shape[0],1))
         
        #diagonal elements should not be zero
        assert np.all(self.nz_index[whichIndices])
        
        #yp is y augmented with zeros           
        yp = np.zeros((self.matrix.shape[0],y.shape[1]), dtype=float, order='C')
        yp[whichIndices, :] = y

        #xp is x augmented with entries that can be ignored
        xp = np.zeros_like(yp)
        A_11 = sliceMatrix(self.matrix, self.nz_index, self.nz_index)
        factor = cholesky(A_11, ordering_method = "natural")

        xp[self.nz_index, :] = factor.solve_Lt(yp[self.nz_index, :], use_LDLt_decomposition = False)
        x = xp[whichIndices, :]
        return x

    