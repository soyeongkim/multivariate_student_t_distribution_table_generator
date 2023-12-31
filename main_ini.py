import numpy as np
from scipy.special import gamma
from scipy.integrate import quad
from scipy.optimize import fsolve
import configparser

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
arr_target_integrity_exp = np.arange(1, 10, step=1)
arr_sample_dof = np.arange(1, 20, step=1)

# init config
config = configparser.ConfigParser()

# set target integrity and num of samples
config.add_section('PL Estimation')
config.set('PL Estimation', 'm_cfg_target_integrity', str(0.0001))
config.set('PL Estimation', 'm_cfg_num_of_samples', str(10))

for integrityIdx in range(len(arr_target_integrity_exp)):
    integrity_dof_table_list = list()
    integrity_k_table_list = list()

    target_integrity = 1 / 10**arr_target_integrity_exp[integrityIdx]
    ini_target_integrity_exp = '1e-' + str(arr_target_integrity_exp[integrityIdx])

    config.add_section(ini_target_integrity_exp)

    for dofIdx in range(len(arr_sample_dof)):
        k = kCalculation(arr_sample_dof[dofIdx], d, target_integrity, 1.0)
        integrity_dof_table_list.append(arr_sample_dof[dofIdx])
        integrity_k_table_list.append(k[0])

    # dof_str = ', '.join(map(str, integrity_dof_table_list))
    # config.set(ini_target_integrity_exp, 'm_cfg_vec_num_of_samples', dof_str)

    k_str = ', '.join(map(str, integrity_k_table_list))
    config.set(ini_target_integrity_exp, 'm_cfg_vec_k', k_str)

    print(integrity_dof_table_list)
    print(integrity_k_table_list)

    with open('multivariate_student_t_distribution_k_dof_table.ini', 'w') as configfile:
        config.write(configfile)