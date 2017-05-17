import numpy as np
import pandas as pd
import sympy as sy
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from math import *
import sys


class Arbitrary:

	def model_format(mode, num_cluster, bound_limit=1e-9):

		data = np.sort(data, axis = 0)
		data_length = len(data)
		raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
		plot_p = np.log(-np.log(1 - raw_probability))

		if mode == 'GM': # Gaussian Mixture
			t = sy.symbols('t')
			u = sy.symarray('u', num_cluster)
			var = sy.symarray('var', num_cluster)
			PDF = [1/(sqrt(2*pi)*j) * sy.exp(((t-i)/(sqrt(2)*j))**2) for i,j in zip(u, var)]

			variable_list = [t, u, var]
			bound = ((bound_limit, 1/bound_limit), (None, None), (None, None))

		elif mode == 'WM': # Weibull Mixture
			t = sy.symbols('t')
			k = sy.symarray('k', num_cluster)
			tor = sy.symarray('tor', num_cluster)
			PDF = sy.diff(1 - sy.exp(- (t/tor)**k), t)
			variable_list = [t, k, tor]
			bound = ((bound_limit, 1/bound_limit), (None, None), (None, None))

		elif mode == 'CM':
			t = sy.symbols('t')
			a = sy.symarray('a', num_cluster)
			b = sy.symarray('b', num_cluster)
			tor = sy.symarray('tor', num_cluster)
			PDF = sy.diff(1 - (1 + 1/a*(t/tor)**(b))**(-a), t)
			variable_list = [t, a, b, tor]
			bound = ((bound_limit, 1/bound_limit), (None, None), (None, None), (bound_limit, 1/bound_limit))

		else:
			print('Mode Error')
			PDF = None
			variable_list = None
			bound = None

		return PDF, variable_list, bound


	def MAP_EST_BFGS(func, var_list, initial_values, bound, tol=1e-9):

		out_values = [[] for k in range(np.size(func[1]))]

		for i in range(np.size(func[1])):
			variable_list = [var_list[x+1][i] for x in range(np.size(var_list,0)-1)]
			object_func = sy.lambdify((variable_list,), func[i], 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=bound)
			out_values[i] = result.x
			func[i] = result.fun

		Total_loglikelihood = sum(func)

		return out_values, Total_loglikelihood


	def EXPA_EST(data, pdf, var_list, initial_values, input_alpha):

		num_cluster = np.size(var_list[1])
		num_data = np.size(data)
		num_para = np.size(var_list,0)

		PDF = pdf[:]

		Q_theta = np.zeros((num_data, num_cluster))

		Q_list = [[] for k in range(num_cluster)] 
		log_likelihood_fun = [[] for k in range(num_cluster)]
		T_likelihood_fun =  [[] for k in range(num_cluster)]
		cluster_list = []
		count = 0

		for x in data:

			for i in range(num_cluster):
				PDF[i] = PDF[i].subs(var_list[0], x)
				for j in range(num_para-1):
					PDF[i] = PDF[i].subs(var_list[j+1][i], initial_values[j][i])

			marginal_prob = np.dot(PDF, input_alpha)
			Q_theta[count][:] = np.array([k * o / marginal_prob for k,o in zip(PDF, input_alpha)])

			for k in range(num_cluster):
				log_likelihood_fun[k].append(sy.log(pdf[k].subs(var_list[0], x) / Q_theta[count][k]) * Q_theta[count][k])
			
			count += 1

		for i in range(num_cluster):
			T_likelihood_fun[i] = sum(log_likelihood_fun[i])

		alpha = [sum(Q_theta[:][i]) / num_data for i in range(num_cluster)]
		
		return T_likelihood_fun, Q_theta, alpha



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








Data_file = pd.ExcelFile(r'C:\Users\Mason\Documents\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709_position.xlsx')  # revise path
#Data_file = pd.ExcelFile(r'C:\Users\Sen\Desktop\1.xlsx')
p_data = Data_file.parse('Sheet1', index_row = None, header = None)
p_data.drop(p_data.columns[[0]], axis = 0, inplace  =True)  # drop first row
p_data = p_data.iloc[:,:].values
data = p_data[:,2]    # data from which column
data = data.astype(np.float32, copy = False)
data = data[~np.isnan(data)]
data = np.array([data])
data = data.T
mean, var, initial_alpha, cluster_list, count = Guassian.EM(data, 2)

print(mean)
print(var)
print(initial_alpha)
print(count)
print(cluster_list)