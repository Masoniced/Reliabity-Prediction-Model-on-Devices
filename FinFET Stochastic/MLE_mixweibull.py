import numpy as np
import pandas as pd
import sympy as sy
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from math import *



def MLE_weibull(data, mode, ratio=1/3, fixa=None, fixb=None, tol=1e-9, custom = None):

	data = np.sort(data, axis = 0)
	data_length = len(data)
	raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
	plot_p = np.log(-np.log(1 - raw_probability))

	if mode=='SC':


		t, k, tor = sy.symbols('t k tor')
		CDF = 1 - sy.exp(- (t/tor)**k)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[7, data[round(data_length*0.63) -1]]
			var_list = [k, tor]
			Probability_bound = ((tol, 1/tol), (min(data),max(data)))
			Likelihood_func = sum([Log_PDF.subs(t, i) for i in data])
		elif fixa==None and fixb!=None:
			initial_values =[1]
			var_list = [k]
			Probability_bound = ((tol,1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (k, fixb)]) for i in data])
		elif fixa!=None and fixb==None:
			initial_values =[data[round(data_length*0.63) -1]]
			var_list = [tor]
			Probability_bound = ((tol,1/tol), (tol,1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (tor, fixa)]) for i in data])


		if custom == None:
			#import pdb; pdb.set_trace()
			object_func = sy.lambdify((var_list,), Likelihood_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = gradient_descent(Likelihood_func, initial, tol)

	elif mode=='BC':

		t, k1, k2, tor1, tor2  = sy.symbols('t k1 k2 tor1 tor2')
		CDF = 1 - sy.exp(- (t/tor1)**k1) * sy.exp(- (t/tor2)**k2)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[ 1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [k1, k2, tor1, tor2]
			Probability_bound = ((tol,1/tol),(tol,1/tol),(tol,1/tol),(tol,1/tol))
			Likelihood_func = sum([Log_PDF.subs(t, i) for i in data])
		else:
			initial_values = [data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [tor1, tor2]
			Probability_bound = ((tol,1/tol), (tol,1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (k1, fixb[0]), (k2, fixb[1])]) for i in data])	 

		if custom == None:
			#import pdb; pdb.set_trace()
			object_func = sy.lambdify((var_list,), Likelihood_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = gradient_descent(Likelihood_func, initial, tol)

	elif mode=='LBC':

		t, k1, k2, tor1, tor2  = sy.symbols('t k1 k2 tor1 tor2')
		CDF = 1 - sy.exp(- (t/tor1)**k1) * sy.exp(- (t/tor2)**k2)

		if fixa==None and fixb==None:
			initial_values =[ 0.001, 0.001, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [k1, k2, tor1, tor2]
			Probability_bound = ((tol,1/tol),(tol,1/tol),(tol,1/tol),(tol,1/tol))
			Least_func = sum([(CDF.subs(t, i) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa==None and fixb!=None:
			initial_values =[1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None))
			Least_func = sum([(CDF.subs([(t, i), (b1, fixb[0]), (b2, fixb[1])]) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa!=None and fixb==None:
			initial_values =[1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [b1, b2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None))
			Least_func = sum([(CDF.subs([(t, i), (a1, fixa[0]), (a2, fixa[1])]) - j)**2 for i,j in zip(data, raw_probability)])
		else:
			initial_values = [data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [tor1, tor2]
			Probability_bound = ((0,None), (0,None))
			Least_func = sum([(CDF.subs([(t, i), (a1, fixa[0]), (a2, fixa[1]), (b1, fixb[0]), (b2, fixb[1])]) - j)**2 for i,j in zip(data, raw_probability)])

		if custom == None:
			object_func = sy.lambdify((var_list,), Least_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = gradient_descent(Least_func, initial, tol)

	elif mode=='LSC':

		t, k, tor = sy.symbols('t k tor')
		CDF = 1 - sy.exp(- (t/tor)**k)


		if fixa==None and fixb==None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [k, tor]
			Probability_bound = ((tol,1/tol),(tol,1/tol))
			Least_func = sum([(CDF.subs(t, i) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa==None and fixb!=None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [a, tor]
			Probability_bound = ((0,None), (0,None))
			Least_func = sum([(CDF.subs[(t, i), (b, fixb)] - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa!=None and fixb==None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [b, tor]
			Probability_bound = ((0,None), (0,None))
			Least_func = sum([(CDF.subs[(t, i), (a, fixa)] - j)**2 for i,j in zip(data, raw_probability)])
		else:
			initial_values =[data[round(data_length*0.63) -1]]
			var_list = [tor]
			Probability_bound = ((0,None))
			Least_func = sum([(CDF.subs[(t, i), (a, fixa), (b, fixb)] - j)**2 for i,j in zip(data, raw_probability)])  

		if custom == None:
			object_func = sy.lambdify((var_list,), Least_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = gradient_descent(Least_func, initial, tol)

	else:
		print("mode error")

	
	if custom == None:

		for i, j in zip(var_list, result.x):
			CDF = CDF.subs(i, j)

		F = sy.lambdify(t, CDF)
		fitting_range = np.linspace(data[0]/1.2, data[-1]*1.2, 200)

		fitting_Prob=[np.log(-np.log(1 - F(i))) for i in fitting_range] 
		return result, data, plot_p, fitting_range, fitting_Prob
	else:
		return MLE_para, error, count, minimum


def MLE(data, mode, ratio=1/3, fixa=None, fixb=None, tol=1e-9, custom = None):

	data = np.sort(data, axis = 0)
	#data = (data-data[0])/(data[-1]-data(0))
	data_length = len(data)
	raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
	plot_p = np.log(-np.log(1 - raw_probability))

	if mode=='SC':


		t, a, b, tor = sy.symbols('t a b tor')
		CDF = 1 - (1 + 1/a*(t/tor)**(b))**(-a)
		Log_PDF = - sy.log(sy.diff(CDF, t))		

		if fixa==None and fixb==None:
			initial_values =[1, 1, data[round(data_length*0.63) -1]]
			var_list = [a, b, tor]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol), (tol, 1/tol))
			Likelihood_func = sum([Log_PDF.subs(t, i) for i in data])
		elif fixa==None and fixb!=None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [a, tor]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (b, fixb)]) for i in data])
		elif fixa!=None and fixb==None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [b, tor]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (a, fixa)]) for i in data])
		else:
			initial_values =[data[round(data_length*0.63) -1]]
			var_list = [tor]
			Probability_bound = ((tol, 1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (a, fixa), (b, fixb)]) for i in data])	

		if custom == None:
			object_func = sy.lambdify((var_list,), Likelihood_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = BFGS(Likelihood_func, initial, tol)

	elif mode=='BC':

		t, a1, a2, b1, b2, tor1, tor2  = sy.symbols('t a1 a2 b1 b2 tor1 tor2')
		CDF = 1 - (1 + 1/a1*(t/tor1)**(b1))**(-a1) * (1 + 1/a2*(t/tor2)**(b2))**(-a2)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[1, 1, 1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, b1, b2, tor1, tor2]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol))
			Likelihood_func = sum([Log_PDF.subs(t, i) for i in data])
		elif fixa==None and fixb!=None:
			initial_values =[1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, tor1, tor2]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (b1, fixb[0]), (b2, fixb[1])]) for i in data])
		elif fixa!=None and fixb==None:
			initial_values =[1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [b1, b2, tor1, tor2]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (a1, fixa[0]), (a2, fixa[1])]) for i in data])
		else:
			initial_values = [data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [tor1, tor2]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (a1, fixa[0]), (a2, fixa[1]), (b1, fixb[0]), (b2, fixb[1])]) for i in data])	 

		if custom == None:
			object_func = sy.lambdify((var_list,), Likelihood_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = BFGS(Likelihood_func, initial, tol)

	elif mode=='LBC':

		t, a1, a2, b1, b2, tor1, tor2 = sy.symbols('t a1 a2 b1 b2 tor1 tor2')
		CDF = 1 - (1 + 1/a1*(t/tor1)**(b1))**(-a1) * (1 + 1/a2*(t/tor2)**(b2))**(-a2)

		if fixa==None and fixb==None:
			initial_values =[1, 1, 1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, b1, b2, tor1, tor2]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol))
			Least_func = sum([(CDF.subs(t, i) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa==None and fixb!=None:
			initial_values =[1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, tor1, tor2]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol))
			Least_func = sum([(CDF.subs([(t, i), (b1, fixb[0]), (b2, fixb[1])]) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa!=None and fixb==None:
			initial_values =[1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [b1, b2, tor1, tor2]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol), (tol, 1/tol), (tol, 1/tol))
			Least_func = sum([(CDF.subs([(t, i), (a1, fixa[0]), (a2, fixa[1])]) - j)**2 for i,j in zip(data, raw_probability)])
		else:
			initial_values = [data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [tor1, tor2]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol))
			Least_func = sum([(CDF.subs([(t, i), (a1, fixa[0]), (a2, fixa[1]), (b1, fixb[0]), (b2, fixb[1])]) - j)**2 for i,j in zip(data, raw_probability)])

		if custom == None:
			object_func = sy.lambdify((var_list,), Least_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = BFGS(Least_func, initial, tol)

	elif mode=='LSC':

		t, a, b, tor = sy.symbols('t a b tor')
		CDF = 1 - (1 + 1/a*(t/tor)**(b))**(-a)

		if fixa==None and fixb==None:
			initial_values =[0.1, 1, data[round(data_length*0.63) -1]]
			var_list = [a, b, tor]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol), (tol, 1/tol))
			Least_func = sum([(CDF.subs(t, i) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa==None and fixb!=None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [a, tor]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol))
			Least_func = sum([(CDF.subs([(t, i), (b, fixb)]) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa!=None and fixb==None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [b, tor]
			Probability_bound = ((tol, 1/tol), (tol, 1/tol))
			Least_func = sum([(CDF.subs([(t, i), (a, fixa)]) - j)**2 for i,j in zip(data, raw_probability)])
		else:
			initial_values =[data[round(data_length*0.63) -1]]
			var_list = [tor]
			Probability_bound = ((tol, 1/tol))
			Least_func = sum([(CDF.subs([(t, i), (a, fixa), (b, fixb)]) - j)**2 for i,j in zip(data, raw_probability)])  

		if custom == None:
			object_func = sy.lambdify((var_list,), Least_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = BFGS(Least_func, initial, tol)

	else:
		print("mode error")


	if custom == None:

		for i, j in zip(var_list, result.x):
			CDF = CDF.subs(i, j)

		F = sy.lambdify(t, CDF)
		#import pdb; pdb.set_trace()
		fitting_range = np.logspace(np.log(data[0]/1.2), np.log(data[-1]*1.2), 200, base=e)


		fitting_Prob=[np.log(-np.log(1 - F(i))) for i in fitting_range]
		
		return result, data, plot_p, fitting_range, fitting_Prob
	
	else:
		return MLE_para, error, count, minimum


def gradient_descent(expr, initial, tol, alpha = 1, mix_iter = 1000):

	var_list = list(initial.keys())
	temp_value = list(initial.values())
	object_func = sy.lambdify((var_list,), expr, 'numpy', dummify=False)	
	gradient = [sy.diff(expr, i) for i in var_list]
	gradient_func = sy.lambdify((var_list,), gradient, 'numpy', dummify=False)

	error = 1
	I_ini = np.identity(len(var_list))
	hessian_inv = I_ini
	count = 0

	while error >= tol:

		g = np.array(gradient_func(temp_value))
		p = - np.dot(hessian_inv, g)
		step = step_checking(object_func, gradient_func, p, temp_value, 1)
		s = step * p
		new_value = temp_value + s
		new_g = np.array(gradient_func(new_value))
		y = new_g - g
		by_product = np.dot(I_ini-np.dot(np.transpose([s]),[y]) / np.dot(s,y), hessian_inv)
		hessian_inv = np.dot(by_product, I_ini-np.dot(np.transpose([y]),[s]) / np.dot(s,y)) + np.dot(np.transpose([s]),[s]) / np.dot(s,y)

		temp_value = new_value
		error = sqrt(np.dot(new_g, new_g))

		count += 1

		if count > mix_iter:
			break


	return dict(zip(var_list, temp_value)), error, count, object_func(temp_value)


def step_checking(f, g, d, x, ini_alpha):

	c1 = 0.01
	c2 = 0.9

	a = 0
	b = inf

	ff = f(x)
	gg = g(x)
	criteria = 1
	alpha = ini_alpha

	while criteria == 1:

		new_x = x + alpha * d
		if any(new_x < 0):
			alpha = alpha / 2
			continue

		new_ff = f(new_x)
		new_gg = g(new_x)

		if new_ff > ff + c1*alpha*np.dot(gg,d):
			b = alpha
			alpha = (a + b) / 2
		elif abs(np.dot(new_gg,d)) > c2*abs(np.dot(gg,d)):
			a = alpha
			if b == inf:
				alpha = 2*a
			else:
				alpha = (a + b) / 2
		else:
			criteria = 0
			
		return alpha




Data_file = pd.ExcelFile(r'C:\Users\Mason\Desktop\GAA.xlsx')  # revise path
#Data_file = pd.ExcelFile(r'C:\Users\Sen\Desktop\1.xlsx')
p_data = Data_file.parse('Sheet1', index_row = None, header = None)
p_data.drop(p_data.columns[[0]], axis = 0, inplace  =True)  # drop first row
p_data = p_data.iloc[:,:].values
data1 = p_data[:,6]    # data from which column
data1 = data1.astype(np.float32, copy = False)
data1 = data1[~np.isnan(data1)]

data2 = p_data[:,7]    # data from which column
data2 = data2.astype(np.float32, copy = False)
data2 = data2[~np.isnan(data2)]

data3 = p_data[:,2]    # data from which column
data3 = data3.astype(np.float32, copy = False)
data3 = data3[~np.isnan(data3)]

data4 = p_data[:,3]    # data from which column
data4 = data4.astype(np.float32, copy = False)
data4 = data4[~np.isnan(data4)]

custom = None
if custom == None:
	# Result1_1, Data1_1, C_Pro1_1, fitting_range1_1, fitting_Prob1_1 = MLE(data1, 'BC', tol=1e-8)
	# Result2_1, Data2_1, C_Pro2_1, fitting_range2_1, fitting_Prob2_1 = MLE_weibull(data1, 'LBC', tol=1e-7)
	# Result3_1, Data3_1, C_Pro3_1, fitting_range3_1, fitting_Prob3_1 = MLE_weibull(data1, 'SC', tol=1e-8)
	# print(Result1_1)
	# print(Result2_1)
	# print(Result3_1)
	# #import pdb; pdb.set_trace()
	# Result1_2, Data1_2, C_Pro1_2, fitting_range1_2, fitting_Prob1_2 = MLE(data2, 'BC', tol=1e-5)
	# Result2_2, Data2_2, C_Pro2_2, fitting_range2_2, fitting_Prob2_2 = MLE_weibull(data2, 'BC', tol=1e-11)
	# Result3_2, Data3_2, C_Pro3_2, fitting_range3_2, fitting_Prob3_2 = MLE_weibull(data2, 'SC', tol=1e-11)
	# print(Result1_2)
	# print(Result2_2)
	# print(Result3_2)

	# Result1_3, Data1_3, C_Pro1_3, fitting_range1_3, fitting_Prob1_3 = MLE(data3, 'BC', tol=1e-11)
	# Result2_3, Data2_3, C_Pro2_3, fitting_range2_3, fitting_Prob2_3 = MLE_weibull(data3, 'BC', tol=1e-11)
	# Result3_3, Data3_3, C_Pro3_3, fitting_range3_3, fitting_Prob3_3 = MLE_weibull(data3, 'SC', tol=1e-11)
	# print(Result1_3)
	# print(Result2_3)
	# print(Result3_3)

	plt.interactive(True)

	# plt.plot(fitting_range2_1, fitting_Prob2_1, 'b-', linewidth=3)
	# plt.plot(fitting_range3_1, fitting_Prob3_1, 'g-', linewidth=3)
	# plt.plot(Data1_1, C_Pro1_1, 'ko', markersize=8)
	# plt.plot(fitting_range1_1, fitting_Prob1_1, 'r-', linewidth=3)

	# plt.plot(fitting_range2_2, fitting_Prob2_2, 'b-', linewidth=3)
	# plt.plot(fitting_range3_2, fitting_Prob3_2, 'g-', linewidth=3)
	# plt.plot(Data1_2, C_Pro1_2, 'ko', markersize=8)
	# plt.plot(fitting_range1_2, fitting_Prob1_2, 'r-', linewidth=3)

	# plt.plot(fitting_range2_3, fitting_Prob2_3, 'b-', linewidth=3)
	# plt.plot(fitting_range3_3, fitting_Prob3_3, 'g-', linewidth=3)
	# plt.plot(Data1_3, C_Pro1_3, 'ko', markersize=8)
	# plt.plot(fitting_range1_3, fitting_Prob1_3, 'r-', linewidth=3)	


	Result1_1, Data1_1, C_Pro1_1, fitting_range1_1, fitting_Prob1_1 = MLE_weibull(data1, 'SC', tol=1e-9)
	Result1_2, Data1_2, C_Pro1_2, fitting_range1_2, fitting_Prob1_2 = MLE_weibull(data2, 'SC', tol=1e-9)
	#import pdb; pdb.set_trace()
	# Result1_3, Data1_3, C_Pro1_3, fitting_range1_3, fitting_Prob1_3 = MLE(data3, 'SC', tol=1e-8)
	# Result1_4, Data1_4, C_Pro1_4, fitting_range1_4, fitting_Prob1_4 = MLE(data4, 'SC', tol=1e-8)
	print(Result1_1)
	print(Result1_2)

	plt.plot(Data1_1, C_Pro1_1, 'ro', markersize=6)
	plt.plot(fitting_range1_1*1.02, fitting_Prob1_1, 'r-', linewidth=3)

	plt.plot(Data1_2, C_Pro1_2, 'go', markersize=6)
	plt.plot(fitting_range1_2, fitting_Prob1_2, 'g-', linewidth=3)

	# plt.plot(Data1_3, C_Pro1_3, 'bo', markersize=6)
	# plt.plot(fitting_range1_3, fitting_Prob1_3, 'b-', linewidth=3)	

	# plt.plot(Data1_4, C_Pro1_4, 'co', markersize=6)
	# plt.plot(fitting_range1_4, fitting_Prob1_4, 'c-', linewidth=3)	

	plt.xscale('log')
	plt.ylim([-4,2])
	plt.xlim([1.5e10,4e10])
	plt.xlabel(r'Time to Failure (s)',{'fontname':'Times New Roman','fontsize':18})
	plt.ylabel(r'ln(-ln(1-F))',{'fontname':'Times New Roman','fontsize':18})

else:
	Result, Error, Iter, Minimum = MLE(data, 'BC', custom = 'self')
	print(Result)
	print(Error)
	print(Iter)
	print(Minimum)