import numpy as np
from scipy.optimize import minimize
import pandas as pd


em = pd.read_csv("salmon-analysis-main/quants_em/sample_2/quant.sf",delimiter="\t")[["Name","TPM"]]
vbem = pd.read_csv("salmon-analysis-main/quants_vbem/sample_2/quant.sf",delimiter="\t")[["Name","TPM"]]


simple_average_results = pd.read_csv("Quant_Joint_Simple_Average.sf",delimiter="\t")["TPM"].to_numpy()


em_data = em["TPM"].to_numpy()
vbem_data = vbem["TPM"].to_numpy()

ground_truth = pd.read_csv("salmon-analysis-main/ground_truth/2.sim.genes.results",delimiter="\t")[["gene_id","TPM"]]
ground_truth = ground_truth.set_index('gene_id').reindex(em['Name']).reset_index()["TPM"].to_numpy()

def MSE(weights,em_iter,vbem_iter,ground_truth):
    weighted_avg = (weights[0] * em_iter + weights[1] * vbem_iter)/np.sum(weights)
    return np.mean((ground_truth - weighted_avg) ** 2)


def MARD(weights,em_iter,vbem_iter,ground_truth):
    weighted_avg = (weights[0] * em_iter + weights[1] * vbem_iter)/np.sum(weights)
    ARD_list = []
    
    for i in range(len(ground_truth)):

        if ground_truth[i] == weighted_avg[i] and ground_truth[i] == 0:
            ARD_list.append(0)
        else:
            ARD_list.append(abs(ground_truth[i] - weighted_avg[i])/(ground_truth[i] + weighted_avg[i]))
    
    return np.mean(ARD_list)


def minimization(intitial_guess, em_iter, vbem_iter, ground_truth, MARD=MARD):
    result = minimize(MARD, x0=intitial_guess, args=(em_iter,vbem_iter,ground_truth),method='L-BFGS-B',bounds=[[0,1],[0,1]])
    return result.x


optimized_weights = minimization([0.5,0.5], em_data, vbem_data, ground_truth, MARD=MARD)
weighted_mean_MARD = MARD(optimized_weights,em_data,vbem_data,ground_truth)


simple_mean_MARD = MARD([0.5,0.5], em_data, vbem_data,ground_truth)



ARD_list = []

for i in range(len(ground_truth)):

    if ground_truth[i] == em_data[i] and ground_truth[i] == 0:
        ARD_list.append(0)
    else:
        ARD_list.append(abs(ground_truth[i] - em_data[i])/(ground_truth[i] + em_data[i]))

em_MARD = np.mean(ARD_list)

ARD_list = []

for i in range(len(ground_truth)):

    if ground_truth[i] == vbem_data[i] and ground_truth[i] == 0:
        ARD_list.append(0)
    else:
        ARD_list.append(abs(ground_truth[i] - vbem_data[i])/(ground_truth[i] + vbem_data[i]))

vbem_MARD = np.mean(ARD_list)


print(f"Weighted Mean MARD: {weighted_mean_MARD}")
print(f"Simple Mean MARD: {simple_mean_MARD}")
print(f"EM MARD: {em_MARD}")
print(f"VBEM MARD: {vbem_MARD}")