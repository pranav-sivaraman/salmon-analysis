import numpy as np
from scipy.optimize import minimize

def MARD(weights,em_iter,vbem_iter):
    aligned_genes = pd.read_csv("salmon-analysis-main/quants_em/sample_1/quant.sf",delimiter="\t")[["Name"]
    
    ground_truth = pd.read_csv("salmon-analysis-main/ground_truth/1.sim.genes.results",delimiter="\t")[["TPM"]].to_numpy()
    ground_truth = ground_truth.set_index('gene_id').reindex(aligned_genes).reset_index()

    weighted_avg = weights[0] * em_iter + weights[1] * vbem_iter
    ARD_list = []
    
    for i in range(len(ground_truth)):

        if ground_truth[i] == weighted_avg[i] and ground_truth[i] == 0:
            ARD_list.append(0)
        else:
            ARD_list.append(abs(ground_truth[i] - weighted_avg[i])/(ground_truth[i] + weighted_avg[i]))
    
    return np.mean(ARD_list)


def minimization(intitial_guess, em_iter, vbem_iter, ground_truth, MARD=MARD):
    result = minimize(MARD, x0=initial_guess, args=(em_iter,vbem_iter,ground_truth), bounds=[(0, 1), (0, 1)])
    return result.x

