'''
author: snowball@USTC
date: 2020.8.31
last update: 2020.8.31
fuction: perform bootstrap analysis for variance
'''

import numpy as np

def linearRegression(x, y): 
    # linear least squares fit 
    # exact same as luminosity.py
    # y = a0 + a1*x
    x = np.array(x); y = np.array(y)
    N = len(x) 
    sumx = sum(x) 
    sumy = sum(y) 
    sumx2 = sum(x**2) 
    sumxy = sum(x*y) 

    # print(sumx); print(sumy)
    # print(N, sumx[0], sumx2[0], sumy[0], sumxy[0])

    try:
        A = np.mat([[N, sumx], [sumx, sumx2]]) 
        b = np.array([sumy, sumxy]) 
    except:
        A = np.mat([[N, float(sumx[0])], [float(sumx[0]), float(sumx2[0])]]) 
        b = np.array([float(sumy[0]), float(sumxy[0])])         

    # A = np.mat([[518, 15530.786], [15530.786, 465840.395]])
    # b = np.array([-782.632, -23494.554])
    a0, a1 = np.linalg.solve(A, b)

    return a0, a1

def bootstrapLF(x, y, m0, m1, times):
    # carry out error estimation through bootstrap method

    temp_m0 = []; temp_m1 = []

    for i in range(times):
        index = np.random.choice(len(x), len(x))
        _x = np.array(x)[index]
        _y = np.array(y)[index]

        _m0, _m1 = linearRegression(_x, _y)
        temp_m0.append((_m0-m0)**2.0); temp_m1.append((_m1-m1)**2.0)

    print('sigma_m0**2 = ' + str(np.sum(temp_m0)/times) + '\n' + 'sigma_m1**2 = ' + \
        str(np.sum(temp_m1)/times))