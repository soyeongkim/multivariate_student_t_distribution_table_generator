import numpy as np
from scipy.special import gamma
from scipy.integrate import quad
from scipy.optimize import fsolve
import configparser

def kCalculation(N, d, target_integrity, init_k):
    # Beta function
    beta_func = gamma((d + N) / 2) / (gamma(d / 2) * gamma(N / 2))

    # Set function for integral
    def tFunc(r):
        return (r**(d - 1)) * (1 + r**2)**(-(N + d) / 2)

    # Object function for finding protection level which satisfies "target integrity"
    # Pillipe Fault paper reference
    def object_function(k):
        integral_result, _ = quad(tFunc, k, np.inf)
        return 2 * beta_func * integral_result - target_integrity

    # Get protection level using "fsolve" function
    options = {'xtol': 1e-20}
    k_solution = fsolve(object_function, init_k, **options)

    return k_solution[0]  # fsolve returns a list, so we take the first element

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
        print(arr_sample_dof[dofIdx], d, arr_target_integrity_exp[integrityIdx],k)
        integrity_k_table_list.append(k)

    k_str = ', '.join(map(str, integrity_k_table_list))
    config.set(ini_target_integrity_exp, 'm_cfg_vec_k', k_str)

    with open('multivariate_student_t_distribution_k_dof_table.ini', 'w') as configfile:
        config.write(configfile)