#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 14:03:45 2021

@author: lucasmurtinho
"""

import numpy as np
import ctypes as ct
from Tree import Node, Tree

import pandas as pd
import time
# import seaborn as sns
# import matplotlib.pyplot as plt

from sklearn.datasets import load_iris, load_digits, load_wine, load_breast_cancer
from sklearn.cluster import KMeans
# from sklearn.datasets import fetch_20newsgroups, fetch_openml, fetch_covtype
# from sklearn.feature_extraction.text import TfidfVectorizer
# from tensorflow.keras.datasets import cifar10

SEEDS = np.arange(1,11)
# SEEDS = [1]
COLUMNS = ['Seed', 'K-means', 'IMM', 'Ex-Greedy']

LIB = ct.CDLL('splitters/best_cut.so')
C_FLOAT_P = ct.POINTER(ct.c_float)
C_INT_P = ct.POINTER(ct.c_int)

LIB.best_cut_single_dim_opt.restype = ct.c_void_p
LIB.best_cut_single_dim_opt.argtypes = [C_FLOAT_P, C_FLOAT_P, C_FLOAT_P, 
                                        C_INT_P, C_INT_P, C_INT_P, 
                                        ct.c_int, ct.c_int, C_FLOAT_P]

def get_distances(data, centers):
    distances = np.zeros((data.shape[0], centers.shape[0]))
    for i in range(centers.shape[0]):
        distances[:,i] = np.linalg.norm(data - centers[i], axis=1) ** 2
    return distances

# def get_possible_cuts(data, centers, dim):
#     possible_cuts = list(data[:,dim])
#     possible_cuts += list(centers[:,dim])
#     return np.unique(possible_cuts)

def get_best_cut_dim(data, valid_data, centers, valid_centers, distances, 
                     dim, func, float_p, int_p):
    n = valid_data.sum()
    k = valid_centers.sum()
    data_order = np.argsort(data[valid_data,dim])
    centers_order = np.argsort(centers[valid_centers,dim])

    data_f = np.asarray(data[valid_data, dim], dtype=np.float64)
    data_p = data_f.ctypes.data_as(float_p)

    centers_f = np.asarray(centers[valid_centers,dim], dtype=np.float64)
    centers_p = centers_f.ctypes.data_as(float_p)

    full_dist_mask = np.outer(valid_data, valid_centers)
    distances_f = np.asarray(distances[full_dist_mask], dtype=np.float64)
    distances_p = distances_f.ctypes.data_as(float_p)
    
    dist_shape = distances_f.reshape(n, k)
    dist_order = np.argsort(dist_shape, axis=1)
    dist_order_f = np.asarray(dist_order, dtype=np.int32).reshape(n*k)
    dist_order_p = dist_order_f.ctypes.data_as(int_p)

    data_order_f = np.asarray(data_order.flatten(), dtype=np.int32)
    data_order_p = data_order_f.ctypes.data_as(int_p)

    centers_order_f = np.asarray(centers_order.flatten(), dtype=np.int32)
    centers_order_p = centers_order_f.ctypes.data_as(int_p)
    
    ans = np.zeros(2, dtype=np.float64)
    ans_p = ans.ctypes.data_as(float_p)
    func(data_p, centers_p, distances_p, dist_order_p, 
         data_order_p, centers_order_p, n, k, ans_p)
    return ans

def best_cut(data, valid_data, centers, valid_centers, distances, 
             verbose=False):
    dim = centers.shape[1]
    best_cut = -np.inf
    best_dim = -1
    best_cost = np.inf
    for i in range(dim):
        # print("for dim", i)
        if len(np.unique(data[valid_data,i])) == 1:
            continue
        ans = get_best_cut_dim(data, valid_data, centers, valid_centers,
                               distances, i, LIB.best_cut_single_dim_opt, 
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

# def test_with_seed(data, k, seed=None):
#     print("for seed {}:".format(seed))
#     start = time.time()
#     km = KMeans(k, random_state=seed)
#     km.fit(data)
#     ks = -km.score(data)
#     end = time.time()
#     print("kmeans done in {:.4f} seconds".format(end-start))
    
#     start = time.time()
#     tree = Tree(k)
#     tree.fit(data, km)
#     ts = tree.score(data)
#     end = time.time()
#     print("IMM done in {:.4f} seconds".format(end-start))
    
#     start = time.time()
#     new_tree = Tree(k)
#     new_tree.tree = fit_tree(data, km.cluster_centers_)
#     ns = new_tree.score(data)
#     end = time.time()
#     print("Ex-Greedy done in {:.4f} seconds".format(end-start))
#     print()
    
#     return ks, ts, ns

# def test_seeds(data, k, seeds):
#     ans = []
#     for s in seeds:
# #         start = time.time()
#         res = [s] + list(test_with_seed(data, k, s))
#         ans.append(res)
# #         end = time.time()
# #         print("done seed {} in {:.4f} seconds".format(s, end-start))
#     return ans

# def build_df(source, columns, data_name):
#     df = pd.DataFrame(source, columns=columns)
#     df['Dataset'] = data_name
#     return df

# digits_data = load_digits().data
# digits_k = 10
# digits_results = test_seeds(digits_data, digits_k, SEEDS)
# df_digits = build_df(digits_results, COLUMNS, 'Digits')
# print(df_digits)