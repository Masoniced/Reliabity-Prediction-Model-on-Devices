import numpy as np
import pandas as pd
from scipy.optimize import minimize
import math
import xlrd


def MLE(data, mode, ratio=1/3, fixa=None, fixb=None, tol=1e-10):

	data = np.sort(data, axis = 0)
	data_length = len(data)
	raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
	plot_p = np.log(-np.log(1 - raw_probability))

	if mode=='SC':

