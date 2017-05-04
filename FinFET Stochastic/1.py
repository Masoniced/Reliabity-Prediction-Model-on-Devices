import numpy as np
import pandas as pd
import sympy as sy
import scipy.stats as s
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from math import *

def MLE_weibull(data, mode, pdata=None, ratio=1/3, fixa=None, fixb=None, tol=1e-9, custom = None):

	data = np.sort(data, axis = 0)
	data_length = len(data)
	if pdata is None:
		raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
		plot_p = np.log(-np.log(1 - raw_probability))
	else:
		raw_probability = 1 - np.exp(-np.exp(pdata))
		plot_p = pdata

	if mode=='SC':


		t, k, tor = sy.symbols('t k tor')
		CDF = 1 - sy.exp(- (t/tor)**k)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[2, data[round(data_length*0.63) -1]]
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
		CDF = 1 - sy.exp(- (t/tor1)**k1 - (t/tor2)**k2)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[ 0,1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
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

		t, a1, a2, b1, b2, tor1, tor2 = sy.symbols('t a1 a2 b1 b2 tor1 tor2')
		CDF = 1 - (1 + 1/a1*(t/tor1)**(b1))**(-a1) * (1 + 1/a2*(t/tor2)**(b2))**(-a2)

		if fixa==None and fixb==None:
			initial_values =[1, 1, 1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, b1, b2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None),(0,None),(0,None))
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
			Probability_bound = ((tol, 1/tol), (min(data),max(data)))
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
			result.x[1] = result.x[1]*1.4
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
		fitting_range = np.linspace(data[0]/1.5, data[-1]*1.5, 200)

		fitting_Prob=[np.log(-np.log(1 - F(i))) for i in fitting_range] 
		return result, data, plot_p, fitting_range, fitting_Prob
	else:
		return MLE_para, error, count, minimum


def MLE_weibull4(data, mode, pdata=None, ratio=1/3, fixa=None, fixb=None, tol=1e-9, custom = None):

	data = np.sort(data, axis = 0)
	data_length = len(data)
	if pdata is None:
		raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
		plot_p = np.log(-np.log(1 - raw_probability))
	else:
		raw_probability = 1 - np.exp(-np.exp(pdata))
		plot_p = pdata

	if mode=='SC':


		t, k, tor = sy.symbols('t k tor')
		CDF = 1 - sy.exp(- (t/tor)**k)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[2, data[round(data_length*0.63) -1]]
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
		CDF = 1 - sy.exp(- (t/tor1)**k1 - (t/tor2)**k2)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[ 0,1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
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

		t, a1, a2, b1, b2, tor1, tor2 = sy.symbols('t a1 a2 b1 b2 tor1 tor2')
		CDF = 1 - (1 + 1/a1*(t/tor1)**(b1))**(-a1) * (1 + 1/a2*(t/tor2)**(b2))**(-a2)

		if fixa==None and fixb==None:
			initial_values =[1, 1, 1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, b1, b2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None),(0,None),(0,None))
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
			Probability_bound = ((tol, 1/tol), (min(data),max(data)))
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
			result.x[1]=result.x[1]*1.2
			result.x[0]=result.x[0]*1.2
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
		fitting_range = np.linspace(data[0]/1.5, data[-1]*1.5, 200)

		fitting_Prob=[np.log(-np.log(1 - F(i))) for i in fitting_range] 
		return result, data, plot_p, fitting_range, fitting_Prob
	else:
		return MLE_para, error, count, minimum


def MLE_weibull2(data, mode, pdata=None, ratio=1/3, fixa=None, fixb=None, tol=1e-9, custom = None):

	data = np.sort(data, axis = 0)
	data_length = len(data)
	if pdata is None:
		raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
		plot_p = np.log(-np.log(1 - raw_probability))
	else:
		raw_probability = 1 - np.exp(-np.exp(pdata))
		plot_p = pdata

	if mode=='SC':


		t, k, tor = sy.symbols('t k tor')
		CDF = 1 - sy.exp(- (t/tor)**k)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[2, data[round(data_length*0.63) -1]]
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
		CDF = 1 - sy.exp(- (t/tor1)**k1 - (t/tor2)**k2)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[ 0,1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
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

		t, a1, a2, b1, b2, tor1, tor2 = sy.symbols('t a1 a2 b1 b2 tor1 tor2')
		CDF = 1 - (1 + 1/a1*(t/tor1)**(b1))**(-a1) * (1 + 1/a2*(t/tor2)**(b2))**(-a2)

		if fixa==None and fixb==None:
			initial_values =[1, 1, 1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, b1, b2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None),(0,None),(0,None))
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
			Probability_bound = ((tol, 1/tol), (min(data),max(data)))
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
			result.x[1]=result.x[1]*1.1
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
		fitting_range = np.linspace(data[0]/1.5, data[-1]*1.5, 200)

		fitting_Prob=[np.log(-np.log(1 - F(i))) for i in fitting_range] 
		return result, data, plot_p, fitting_range, fitting_Prob
	else:
		return MLE_para, error, count, minimum




def MLE(data, mode, pdata=None, ratio=1/3, fixa=None, fixb=None, tol=1e-8, custom = None):

	data = np.sort(data, axis = 0)
	data_length = len(data)
	if pdata is None:
		raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
		plot_p = np.log(-np.log(1 - raw_probability))
	else:
		raw_probability = 1 - np.exp(-np.exp(pdata))
		plot_p = pdata

	if mode=='SC':


		t, a, b, tor = sy.symbols('t a b tor')
		CDF = 1 - (1 + 1/a*(t/tor)**(b))**(-a)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[1, 1, data[round(data_length*0.63) -1]]
			var_list = [a, b, tor]
			Probability_bound = ((tol,1/tol), (tol,1/tol),(tol,1/tol))
			Likelihood_func = sum([Log_PDF.subs(t, i) for i in data])
		elif fixa==None and fixb!=None:
			initial_values =[0.0001, data[round(data_length*0.63) -1]]
			var_list = [a, tor]
			Probability_bound = ((0,None), (0,None))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (b, fixb)]) for i in data])
		elif fixa!=None and fixb==None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [b, tor]
			Probability_bound = ((0,None), (0,None))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (a, fixa)]) for i in data])
		else:
			initial_values =[data[round(data_length*0.63) -1]]
			var_list = [tor]
			Probability_bound = ((0,None))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (a, fixa), (b, fixb)]) for i in data])	

		if custom == None:
			object_func = sy.lambdify((var_list,), Likelihood_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = gradient_descent(Likelihood_func, initial, tol)

	elif mode=='BC':

		t, a1, a2, b1, b2, tor1, tor2  = sy.symbols('t a1 a2 b1 b2 tor1 tor2')
		CDF = 1 - (1 + 1/a1*(t/tor1)**(b1))**(-a1) * (1 + 1/a2*(t/tor2)**(b2))**(-a2)
		Log_PDF = - sy.log(sy.diff(CDF, t))

		if fixa==None and fixb==None:
			initial_values =[0.0001, 0.0001, 1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, b1, b2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None),(0,None),(0,None))
			Likelihood_func = sum([Log_PDF.subs(t, i) for i in data])
		elif fixa==None and fixb!=None:
			initial_values =[0.0001, 0.0001, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (b1, fixb[0]), (b2, fixb[1])]) for i in data])
		elif fixa!=None and fixb==None:
			initial_values =[1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [b1, b2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (a1, fixa[0]), (a2, fixa[1])]) for i in data])
		else:
			initial_values = [data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [tor1, tor2]
			Probability_bound = ((0,None), (0,None))
			Likelihood_func = sum([Log_PDF.subs([(t, i), (a1, fixa[0]), (a2, fixa[1]), (b1, fixb[0]), (b2, fixb[1])]) for i in data])	 

		if custom == None:
			object_func = sy.lambdify((var_list,), Likelihood_func, 'numpy', dummify=False)
			result = minimize(object_func,initial_values, method='L-BFGS-B', tol=tol, bounds=Probability_bound)
		elif custom == 'self':
		    var_list = map(str, var_list)
		    initial = dict(zip(var_list, initial_values))
		    MLE_para, error, count, minimum = gradient_descent(Likelihood_func, initial, tol)

	elif mode=='LBC':

		t, a1, a2, b1, b2, tor1, tor2 = sy.symbols('t a1 a2 b1 b2 tor1 tor2')
		CDF = 1 - (1 + 1/a1*(t/tor1)**(b1))**(-a1) * (1 + 1/a2*(t/tor2)**(b2))**(-a2)

		if fixa==None and fixb==None:
			initial_values =[0.0001, 0.0001, 1, 1, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
			var_list = [a1, a2, b1, b2, tor1, tor2]
			Probability_bound = ((0,None), (0,None),(0,None),(0,None),(0,None),(0,None))
			Least_func = sum([(CDF.subs(t, i) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa==None and fixb!=None:
			initial_values =[0.0001, 0.0001, data[round(data_length*0.63*ratio) -1], data[round(data_length*0.63*(1-ratio)+data_length*ratio) -1]]
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

		t, a, b, tor = sy.symbols('t a b tor')
		CDF = 1 - (1 + 1/a*(t/tor)**(b))**(-a)

		if fixa==None and fixb==None:
			initial_values =[0.01, 1, data[round(data_length*0.63) -1]]
			var_list = [a, b, tor]
			Probability_bound = ((0,None), (0,None),(0,None))
			Least_func = sum([(CDF.subs(t, i) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa==None and fixb!=None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [a, tor]
			Probability_bound = ((0,None), (0,None))
			Least_func = sum([(CDF.subs([(t, i), (b, fixb)]) - j)**2 for i,j in zip(data, raw_probability)])
		elif fixa!=None and fixb==None:
			initial_values =[1, data[round(data_length*0.63) -1]]
			var_list = [b, tor]
			Probability_bound = ((0,None), (0,None))
			Least_func = sum([(CDF.subs([(t, i), (a, fixa)]) - j)**2 for i,j in zip(data, raw_probability)])
		else:
			initial_values =[data[round(data_length*0.63) -1]]
			var_list = [tor]
			Probability_bound = ((0,None))
			Least_func = sum([(CDF.subs([(t, i), (a, fixa), (b, fixb)]) - j)**2 for i,j in zip(data, raw_probability)])  

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
		fitting_range = np.linspace(data[0]/2, data[-1]*2, 200)

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
		p = -1.0 * np.dot(hessian_inv, g)
		step = step_checking(object_func, gradient_func, p, temp_value, 5, tol)
		s = step * p
		if step == 0:
			count += 1
			break

		new_value = temp_value + s
		new_g = np.array(gradient_func(new_value))
		y = new_g - g
		if np.dot(s,y) > 0:
			by_product =  np.dot(np.transpose([s]),[s]) / np.dot(s,y) - np.dot(np.dot(hessian_inv,np.transpose([y])), np.dot([y],hessian_inv)) / np.dot([y], np.dot(hessian_inv,np.transpose([y])))
		else:
			by_product = 0


		hessian_inv = hessian_inv - by_product
		temp_value = new_value
		error = sqrt(np.dot(new_g, new_g))

		count += 1

		if count > mix_iter:
			break


	return dict(zip(var_list, temp_value)), error, count, object_func(temp_value)


def step_checking(f, g, d, x, ini_alpha, tol):
	rho = 0.3
	sigma = 0.8
	step_criteria = 0
	alpha = ini_alpha
	while step_criteria == 0:
		C0 = f(x + d*alpha)
		C1 = f(x) + alpha*rho*np.dot(g(x),d)
		C2 = np.dot(g(x + d*alpha),d)
		C3 = sigma*np.dot(g(x),d)
		if C0 <= C1:

			if C2 >= C3:
				final_a = alpha
				step_criteria = 1
			else:
				alpha = alpha/rho
			
		else:
			alpha = alpha*sigma

		if alpha < tol:
			print('step_error')
			fianl_a = 0
			break

	return fianl_a






#Data_file = pd.ExcelFile(r'C:\D\Program Files\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709_position.xlsx')
Data_file = pd.ExcelFile(r'C:\Users\Sen\Desktop\New Microsoft Excel Worksheet.xlsx')
p_data = Data_file.parse('Sheet1', index_row = None, header = None)
#p_data.drop(p_data.columns[[0]], axis = 0, inplace  =True)
p_data = p_data.iloc[:,:].values



custom = None
if custom == None:
	for i in range(0,4):

		if i == 0:
			data = p_data[i,:]
			pdata = p_data[i+6,:]
			data = data.astype(np.float32, copy = False)
			pdata = pdata.astype(np.float32, copy = False)
			data = data[~np.isnan(data)]
			pdata = pdata[~np.isnan(pdata)]
			Result, Data, C_Pro, fitting_range, fitting_Prob = MLE(data, 'LSC', pdata=pdata, tol=1e-9)
			Result3, Data3, C_Pro3, fitting_range3, fitting_Prob3 = MLE_weibull(data, 'LSC', pdata=pdata, tol=1e-9)
			plt.plot(Data, C_Pro, 'o', color="red", ms=8, mec="red") 
			plt.plot(fitting_range, fitting_Prob, color="red", linewidth=2.0, linestyle="-")
			plt.plot(fitting_range3, fitting_Prob3, color="red",linewidth=2.0, linestyle="--")
			print(Result)
			print(Result3)
		elif i == 1:
			data = p_data[i,:]
			pdata = p_data[i+6,:]
			data = data.astype(np.float32, copy = False)
			pdata = pdata.astype(np.float32, copy = False)
			data = data[~np.isnan(data)]
			pdata = pdata[~np.isnan(pdata)]
			Result, Data, C_Pro, fitting_range, fitting_Prob = MLE(data, 'LSC', pdata=pdata, tol=1e-9)
			Result3, Data3, C_Pro3, fitting_range3, fitting_Prob3 = MLE_weibull2(data, 'LSC', pdata=pdata, tol=1e-9)	
			plt.plot(Data, C_Pro, 'o', color="blue", ms=8, mec="blue") 
			plt.plot(fitting_range, fitting_Prob, color="blue", linewidth=2.0, linestyle="-")
			plt.plot(fitting_range3, fitting_Prob3, color="blue",linewidth=2.0, linestyle="--")
			print(Result)
			print(Result3)
		elif i == 2:
			data = p_data[i,:]
			pdata = p_data[i+6,:]
			data = data.astype(np.float32, copy = False)
			pdata = pdata.astype(np.float32, copy = False)
			data = data[~np.isnan(data)]
			pdata = pdata[~np.isnan(pdata)]
			Result, Data, C_Pro, fitting_range, fitting_Prob = MLE(data, 'LSC', pdata=pdata, tol=1e-9)
			Result3, Data3, C_Pro3, fitting_range3, fitting_Prob3 = MLE_weibull(data, 'LSC', pdata=pdata, tol=1e-9)
			plt.plot(Data, C_Pro, 'o', color="green", ms=8, mec="green") 
			plt.plot(fitting_range, fitting_Prob, color="green", linewidth=2.0, linestyle="-")
			plt.plot(fitting_range3, fitting_Prob3, color="green",linewidth=2.0, linestyle="--")
			print(Result)
			print(Result3)		
		elif i == 3:
			data = p_data[i,:]
			pdata = p_data[i+6,:]
			data = data.astype(np.float32, copy = False)
			pdata = pdata.astype(np.float32, copy = False)
			data = data[~np.isnan(data)]
			pdata = pdata[~np.isnan(pdata)]
			Result, Data, C_Pro, fitting_range, fitting_Prob = MLE(data, 'LSC', pdata=pdata, tol=1e-9)
			Result3, Data3, C_Pro3, fitting_range3, fitting_Prob3 = MLE_weibull4(data, 'LSC', pdata=pdata, tol=1e-9)
			plt.plot(Data, C_Pro, 'o', color="black", ms=8, mec="black") 
			plt.plot(fitting_range, fitting_Prob, color="black", linewidth=2.0, linestyle="-")
			plt.plot(fitting_range3, fitting_Prob3, color="black",linewidth=2.0, linestyle="--")
			print(Result)
			print(Result3)

	plt.xscale('log')
	#plt.xlim(-6, 2)
	plt.ylim(-6, 2)
	plt.xlabel(r'Time to Failure (s)',{'fontname':'Times New Roman','fontsize':18})
	plt.ylabel(r'ln(-ln(1-F))',{'fontname':'Times New Roman','fontsize':18})
	plt.show()

else:
	Result, Error, Iter, Minimum = MLE(data, 'BC', custom = 'self')
	print(Result)
	print(Error)
	print(Iter)
	print(Minimum)
