import numpy as np
from scipy.special import gamma
from scipy.integrate import quad
from scipy.optimize import fsolve
import csv

def kCalculation(N, d, target_integrity, init_k):
    result_k = fsolve(objFunc, init_k, args=(N, d, target_integrity), xtol=1e-1000, maxfev=10000)

    return result_k

def objFunc(k, N, d, target_integrity):
    return 2 * betaFunc(N,d) * quadInteger(k, N, d) - target_integrity

def tFunc(r, N, d):
    return (pow(r,(d-1))) * pow((1 + pow(r,2)),(-(N+d)/2))

def betaFunc(d, N):
    return gamma((d+N)/2) / (gamma(d/2)*gamma(N/2))

def quadInteger(k, N, d):
    y, err = quad(tFunc, k, np.infty, args=(N, d))
    return y

# fixed param
d = 2.0 # dimension

# non-fixed param
arr_target_integrity = np.array([0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001, 0.00000000001])
arr_sample_dof = np.arange(1, 20, step=1) # 1:1:20

# Calculate 'k' value using solver
whole_integrity_dof_k_table = list()

for integrityIdx in range(len(arr_target_integrity)):
    integrity_dof_k_table_list = list()

    # save in csv file
    file_name_tir = str(arr_target_integrity[integrityIdx]) + '.csv'
    file_name = open(file_name_tir, "w", newline="")
    writer = csv.writer(file_name)

    for dofIdx in range(len(arr_sample_dof)):
        k = kCalculation(arr_sample_dof[dofIdx], d, arr_target_integrity[integrityIdx], 1.0)
        integrity_dof_k_table_list.append([arr_sample_dof[dofIdx], k[0]])
    
    print(integrity_dof_k_table_list)
    writer.writerows(integrity_dof_k_table_list)
    file_name.close()