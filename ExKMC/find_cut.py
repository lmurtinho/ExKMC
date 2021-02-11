#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 20:17:38 2021

@author: lucasmurtinho
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:03:45 2021

@author: lucasmurtinho
"""

import numpy as np
import ctypes as ct
from Tree import Node

SEEDS = np.arange(1,11)
COLUMNS = ['Seed', 'K-means', 'IMM', 'Ex-Greedy']

LIB = ct.CDLL('splitters/best_cut.so')
C_FLOAT_P = ct.POINTER(ct.c_float)
C_INT_P = ct.POINTER(ct.c_int)

LIB.best_cut_single_dim.restype = ct.c_void_p
LIB.best_cut_single_dim.argtypes = [C_FLOAT_P, C_FLOAT_P, C_FLOAT_P, 
                                    C_INT_P, ct.c_int, ct.c_int, 
                                    C_FLOAT_P]

def get_distances(data, centers):
    distances = np.zeros((data.shape[0], centers.shape[0]))
    for i in range(centers.shape[0]):
        distances[:,i] = np.linalg.norm(data - centers[i], axis=1) ** 2
    return distances

def get_best_cut_dim(data, valid_data, centers, valid_centers, 
                     distances_pointer, dist_order_pointer,
                     n, k, dim, func, float_p, int_p):

    data_f = np.asarray(data[valid_data, dim], dtype=np.float64)
    data_p = data_f.ctypes.data_as(float_p)

    centers_f = np.asarray(centers[valid_centers,dim], dtype=np.float64)
    centers_p = centers_f.ctypes.data_as(float_p)
    
    ans = np.zeros(2, dtype=np.float64)
    ans_p = ans.ctypes.data_as(float_p)
    func(data_p, centers_p, distances_pointer, 
         dist_order_pointer, n, k, ans_p)
    return ans

def best_cut(data, valid_data, centers, valid_centers, distances, 
             verbose=False):
    dim = centers.shape[1]
    best_cut = -np.inf
    best_dim = -1
    best_cost = np.inf
    
    n = valid_data.sum()
    k = valid_centers.sum()
    
    full_dist_mask = np.outer(valid_data, valid_centers)
    distances_f = np.asarray(distances[full_dist_mask], dtype=np.float64)
    distances_p = distances_f.ctypes.data_as(C_FLOAT_P)
    
    dist_shape = distances_f.reshape(n, k)
    dist_order = np.argsort(dist_shape, axis=1)
    dist_order_f = np.asarray(dist_order, dtype=np.int32).reshape(n*k)
    dist_order_p = dist_order_f.ctypes.data_as(C_INT_P)
    
    for i in range(dim):
        if len(np.unique(data[valid_data,i])) == 1:
            continue
        ans = get_best_cut_dim(data, valid_data, centers, valid_centers,
                               distances_p, dist_order_p, n, k, i, 
                               LIB.best_cut_single_dim, 
                               C_FLOAT_P, C_INT_P)
        cut, cost = ans
        if verbose:
            print("for dimension {}: {}, {}".format(i, cut, cost))
        if cost < best_cost:
            best_cut = cut
            best_dim = i
            best_cost = cost
        if verbose:
            print()
    return best_dim, best_cut, best_cost

def build_tree(data, centers, distances, valid_centers, valid_data):
    node = Node()
    if valid_centers.sum() == 1:
        node.value = np.argmax(valid_centers)
        return node
    dim, cut, cost = best_cut(data, valid_data, centers, 
                              valid_centers, distances)
    node.feature = dim
    node.value = cut
    
    n = data.shape[0]
    left_valid_data = np.zeros(n, dtype=bool)
    right_valid_data = np.zeros(n, dtype=bool)
    for i in range(n):
        if valid_data[i]:
            if data[i,dim] <= cut:
                left_valid_data[i] = True
            else:
                right_valid_data[i] = True

    k = centers.shape[0]
    left_valid_centers = np.zeros(k, dtype=bool)
    right_valid_centers = np.zeros(k, dtype=bool)
    for i in range(k):
        if valid_centers[i]:
            if centers[i, dim] <= cut:
                left_valid_centers[i] = True
            else:
                right_valid_centers[i] = True
    
    node.left = build_tree(data, centers, distances, 
                           left_valid_centers, left_valid_data)
    node.right = build_tree(data, centers, distances, 
                            right_valid_centers, right_valid_data)
    
    return node

def fit_tree(data, centers):
    k, d = centers.shape
    n = data.shape[0]
    valid_centers = np.ones(k, dtype=bool)
    valid_data = np.ones(n, dtype=bool)
    distances = get_distances(data, centers)
    return build_tree(data, centers, distances, valid_centers, valid_data)