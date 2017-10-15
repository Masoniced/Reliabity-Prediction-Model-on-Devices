import numpy as np
import pandas as pd
import sympy as sy
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from math import *
import sys


class Guassian:

	#data, initial_mean, initial_var as matrix
	def EM(data, num_cluster, initial_mean=None, initial_var=None, tol=1e-9):

		num_data, dim_data = np.shape(data)

		if initial_mean == None:

			mean = np.zeros((num_cluster, dim_data))
			ini_v = np.zeros((1,dim_data))
			for i in range(dim_data):
				ini_v[i] = np.std(data[:,i])
				for j in range(num_cluster):
					ind = np.around(np.random.uniform(0,1)*num_data)
					mean[j,i] = np.mean(data[ind.astype(int):,i])

		if initial_var == None:
			
			var = np.zeros((dim_data, dim_data, num_cluster))
			for i in range(num_cluster):
				cluster_v = np.eye(dim_data)
				np.fill_diagonal(cluster_v,ini_v)
				var[:,:,i] = cluster_v

		criteria = 1
		initial_alpha = np.ones((num_cluster,1)) / num_cluster
		Q_theta = np.zeros((num_data, num_cluster))
		Total_likelihood = 0
		Old_total = 1
		cluster_list = np.zeros((num_data,1))
		count = 0


		while criteria == 1:

			count += 1

			# E-step
			for i in range(num_data):

				Q_theta[[i],:] = np.array([exp(-(data[i,:] - mean[k,:]) * np.linalg.inv(var[:,:,k]) * np.transpose([data[i,:] - mean[k,:]]) / 2) / ((2*pi)**(num_cluster/2) * (np.linalg.det(var[:,:,k]))**(1/2)) * initial_alpha[k] for k in range(num_cluster)]).T
				cluster_total = sum(Q_theta[i,:])
				Total_likelihood += cluster_total
				if cluster_total == 0:
					cluster_total = 1
				Q_theta[i,:] = Q_theta[i,:] / cluster_total

			# Check convergence
			#import pdb; pdb.set_trace()
			if Total_likelihood == 0:
				criteria = 0
				break
			if abs((Old_total - Total_likelihood)/Old_total) < tol:
				criteria = 0
				break

			# M-step
			initial_alpha = np.array([sum(Q_theta[:,k]) / num_data for k in range(num_cluster)]) # Alpha

			for i in range(num_cluster):
				mean[i,:] = np.array([np.dot(Q_theta[:,i], data[:,k]) / (initial_alpha[i]*num_data) for k in range(dim_data)]) # mean

			for i in range(num_cluster):
				var[:,:,i] = sum([np.transpose([data[k,:] - mean[i,:]]) * (data[k,:] - mean[i,:]) * Q_theta[k,i] for k in range(num_data)]) / (initial_alpha[i]*num_data)  # variance

			Old_total = Total_likelihood
			Total_likelihood = 0

		for i in range(num_data):
			mv, ind = max([(mv, ind) for ind, mv in enumerate(Q_theta[i,:])])
			cluster_list[i] = ind


		return mean, var, initial_alpha, cluster_list, count








Data_file = pd.ExcelFile(r'C:\Users\Mason\Documents\Project')  # revise path
p_data = Data_file.parse('Sheet1', index_row = None, header = None)
p_data.drop(p_data.columns[[0]], axis = 0, inplace  =True)  # drop first row
p_data = p_data.iloc[:,:].values
data = p_data[:,2]    # data from which column
data = data.astype(np.float32, copy = False)
data = data[~np.isnan(data)]
data = np.array([data])
data = data.T
mean, var, initial_alpha, cluster_list, count = Guassian.EM(data, 4)

print(mean)
print(var)
print(initial_alpha)
print(count)
print(cluster_list)