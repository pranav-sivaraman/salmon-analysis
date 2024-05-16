import numpy as np
import pandas as pd
import tqdm
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error


def MARD(weights,em_iter,vbem_iter,ground_truth):
    weighted_avg = (weights[0] * em_iter + weights[1] * vbem_iter)
    ARD_list = []
    for i in range(len(ground_truth)):

        if ground_truth[i] == weighted_avg[i] and ground_truth[i] == 0:
            ARD_list.append(0)
        else:
            ARD_list.append(np.abs(ground_truth[i] - weighted_avg[i])/(ground_truth[i] + weighted_avg[i]))
    
    return np.mean(ARD_list)



em = pd.read_csv("salmon-analysis-main/quants_em/sample_1/quant.sf",delimiter="\t")[["Name","TPM"]]
# vbem = pd.read_csv("salmon-analysis-main/quants_vbem/sample_20/quant.sf",delimiter="\t")[["Name","TPM"]]


em_data = np.fromfile("intermediate_alphas_em.bin",dtype=np.double)
vbem_data = np.fromfile("intermediate_alphas_vbem.bin",dtype=np.double)

em_data = em_data.reshape((1000,int(em_data.shape[0]/1000)))
vbem_data = vbem_data.reshape((1000,int(vbem_data.shape[0]/1000)))

ground_truth = pd.read_csv("salmon-analysis-main/ground_truth/1.sim.genes.results",delimiter="\t")[["gene_id","TPM"]]
ground_truth = ground_truth.set_index('gene_id').reindex(em['Name']).reset_index()



# def mse(weights): 
#     weighted_avg = (weights[0] * em_iter + weights[1] * vbem_iter)
#     return np.mean((weighted_avg - ground_truth["TPM"].to_numpy()[0:1000])**2)


initial_guess = [0.1,0.1]
min_iter = min(em_data.shape[1], vbem_data.shape[1])
weight_np = None
mse_np = np.array([])

for i in range(em_data.shape[1]):
    # for j in range(vbem_data):
    em_iter = em_data[:,i]/np.sum(em_data[:,i])
    vbem_iter = em_data[:,i]/np.sum(vbem_data[:,i])

    result = minimize(MARD, initial_guess, bounds=[(0, 1), (0, 1)])
    

    weighted_avg = (result.x[0] * em_iter + result.x[1] * vbem_iter)
    mse_val = np.mean((weighted_avg - ground_truth["TPM"].to_numpy()[0:1000])**2)

    mse_np = np.append(mse_np, mse_val)

    if weight_np is None:
        weight_np = np.array([result.x])
    
    else:
        weight_np = np.concatenate((weight_np, [result.x]))



# result = minimize(mse, initial_guess, bounds=[(0, 1), (0, 1)])

# optimized_weights = result.x
# weighted_avg = (optimized_weights[0] * em["TPM"].to_numpy() + optimized_weights[1] * vbem["TPM"].to_numpy())/np.sum(optimized_weights)

# mse_val = np.mean((weighted_avg - ground_truth["TPM"].to_numpy())**2)

# # plt.plot(em["TPM"][1:25],vbem["TPM"][1:25])

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Plot the 3D scatter plot
# ax.scatter(em["TPM"][1:25], vbem["TPM"][1:25], ground_truth["TPM"][1:25])
# ax.set_xlabel("EM")
# ax.set_ylabel("VBEM")
# ax.set_zlabel("Ground Truth")

# plt.show()


