#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:26:56 2021

@author: lucasmurtinho
"""

import numpy as np
import ctypes as ct
from Tree import Node

LIB = ct.CDLL('splitters/best_cut.so')
C_FLOAT_P = ct.POINTER(ct.c_float)
C_INT_P = ct.POINTER(ct.c_int)

LIB.best_cut_single_dim.restype = ct.c_void_p
LIB.best_cut_single_dim.argtypes = [C_FLOAT_P, C_FLOAT_P, C_FLOAT_P, C_INT_P,
                                    C_FLOAT_P, C_FLOAT_P, C_INT_P, C_INT_P, 
                                    ct.c_int, ct.c_int, ct.c_int, C_FLOAT_P]

def get_distances(data, centers):
    distances = np.zeros((data.shape[0], centers.shape[0]))
    for i in range(centers.shape[0]):
        distances[:,i] = np.linalg.norm(data - centers[i], axis=1) ** 2
    return distances

def get_possible_cuts(data, centers, dim):
    min_cut = centers[:,dim].min()
    max_cut = centers[:,dim].max()
    possible_cuts = list(data[:,dim])
    possible_cuts += list(centers[:,dim])
    # CHECK: no need for filter below because it's done in best_cut function
    possible_cuts = [i for i in possible_cuts 
                     if (i >= min_cut) and (i < max_cut)]
    return np.unique(possible_cuts)

def get_best_cut_dim(data, valid_data, centers, valid_centers, distances, 
                     dim, cuts, func, float_p, int_p):
    # print(valid_data.sum(), "valid data")
    n = valid_data.sum()
    k = valid_centers.sum()
    # costs = distances[valid_data].min(axis=1)
    data_order = np.argsort(data[valid_data,dim])
    centers_order = np.argsort(centers[valid_centers,dim])
    n_cuts = cuts.shape[0]

    data_f = np.asarray(data[valid_data, dim], dtype=np.float64)
    # data_f = data.flatten()
    data_p = data_f.ctypes.data_as(float_p)

    centers_f = np.asarray(centers[valid_centers,dim], dtype=np.float64)
    # centers_f = centers.flatten()
    centers_p = centers_f.ctypes.data_as(float_p)

    full_dist_mask = np.outer(valid_data, valid_centers)
    distances_f = np.asarray(distances[full_dist_mask], dtype=np.float64)
    # distances_f = distances.flatten()
    distances_p = distances_f.ctypes.data_as(float_p)
    
    dist_shape = distances_f.reshape(n, k)
    dist_order = np.argsort(dist_shape, axis=1)
    dist_order_f = np.asarray(dist_order, dtype=np.int32)
    # dist_order_f = dist_order[valid_data].flatten()
    dist_order_p = dist_order_f.ctypes.data_as(int_p)

    cuts_f = np.asarray(cuts.flatten(), dtype=np.float64)
    # cuts_f = cuts.flatten()
    cuts_p = cuts_f.ctypes.data_as(float_p)

    costs_f = np.asarray(distances_f.reshape(n,k).min(axis=1), 
                         dtype=np.float64)
    # costs_f = np.asarray(costs.flatten(), dtype=np.float64)
    # costs_f = costs.flatten()
    costs_p = costs_f.ctypes.data_as(float_p)

    data_order_f = np.asarray(data_order.flatten(), dtype=np.int32)
    # data_order_f = data_order.flatten()
    data_order_p = data_order_f.ctypes.data_as(int_p)

    centers_order_f = np.asarray(centers_order.flatten(), dtype=np.int32)
    # centers_order_f = centers_order.flatten()
    centers_order_p = centers_order_f.ctypes.data_as(int_p)
    
    ans = np.zeros(2, dtype=np.float64)
    ans_p = ans.ctypes.data_as(float_p)
    # print(data_f.shape, centers_f.shape, distances_f.shape, 
    #       dist_order_f.shape, cuts_f.shape, costs_f.shape, data_order_f.shape,
    #       centers_order_f.shape, ans.shape)
    # print("enter C")
    func(data_p, centers_p, distances_p, dist_order_p, cuts_p, costs_p,
         data_order_p, centers_order_p, n, k, n_cuts, ans_p)
    # print("left C")
    return ans

def best_cut(data, valid_data, centers, valid_centers, distances, 
             possible_cuts, verbose=False):
    dim = centers.shape[1]
    best_cut = -np.inf
    best_dim = -1
    best_cost = np.inf
    for i in range(dim):
        if len(np.unique(data[valid_data,i])) == 1:
            continue
        cuts = possible_cuts[i].copy()
        min_cut = centers[valid_centers,i].min()
        max_cut = centers[valid_centers,i].max()
        cuts = cuts[cuts >= min_cut]
        cuts = cuts[cuts < max_cut]
        # print(len(cuts), "cuts")
        if len(cuts):
            ans = get_best_cut_dim(data, valid_data, centers, valid_centers,
                                   distances, i, cuts, LIB.best_cut_single_dim, 
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

def build_tree(data, centers, possible_cuts, distances, valid_centers, 
               valid_data):
    # print("{} valid centers".format(valid_centers.sum()))
    node = Node()
    if valid_centers.sum() == 1:
        node.value = np.argmax(valid_centers)
        return node
    dim, cut, cost = best_cut(data, valid_data, centers, valid_centers, 
                              distances, possible_cuts)
    # print(dim, cut, cost)
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
    node.left = build_tree(data, centers, possible_cuts, distances, 
                           left_valid_centers, left_valid_data)
    node.right = build_tree(data, centers, possible_cuts, distances, 
                            right_valid_centers, right_valid_data)
    
    return node

def fit_tree(data, centers):
    k, d = centers.shape
    n = data.shape[0]
    valid_centers = np.ones(k, dtype=bool)
    valid_data = np.ones(n, dtype=bool)
    possible_cuts = [get_possible_cuts(data, centers, i) 
                     for i in range(d)]
    distances = get_distances(data, centers)
    return build_tree(data, centers, possible_cuts, distances, valid_centers, 
                      valid_data)