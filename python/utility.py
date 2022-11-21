#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 10:46:43 2022

@author: psalehin
"""

def sliceMatrix(A, l1, l2):
    A1 = A[l1,:]
    A12= A1[:,l2]
    return A12