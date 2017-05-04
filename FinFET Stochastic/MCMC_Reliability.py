import numpy as np
import pandas as pd
import sympy as sy
from scipy.optimize import minimize
from math import *
import xlrd

def sampler(data, samples=4, mu_init=.5, proposal_width=.5, plot=False, mu_prior_mu=0, mu_prior_sd=1.):
    mu_current = mu_init
    posterior = [mu_current]
    for i in range(samples):
        # suggest new position
        mu_proposal = norm(mu_current, proposal_width).rvs()

        # Compute likelihood by multiplying probabilities of each data point
        likelihood_current = norm(mu_current, 1).pdf(data).prod()
        likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
        
        # Compute prior probability of current and proposed mu        
        prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
        posterior_proposal =  norm_hist(mu_prososal)
        prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
        
        p_current = likelihood_current * prior_current
        p_proposal = likelihood_proposal * prior_proposal
        
        # Accept proposal?
        p_accept = p_proposal / p_current
        
        # Usually would include prior probability, which we neglect here for simplicity
        accept = np.random.rand() < p_accept
        
        if plot:
            plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accept, posterior, i)
        
        if accept:
            # Update position
            mu_current = mu_proposal
        
        posterior.append(mu_current)
        
    return posterior




file_name = r"C:\Users\Sen\Desktop\Raw-Curves Files\Logan-7C_w12_TDDB_25C_Compiled Raw.txt"
File = pd.read_csv(file_name, sep ='\t', header = 0)
Columns = File.columns
Result = pd.DataFrame(columns = Columns)
Ini_key = File.iat[0,2]
criteria = 0
Low_resistance = 1000

for i in range(1, len(File.index)):
    if File.iat[i,2] != Ini_key:
        Ini_key = File.iat[i,2]
        criteria = 0
    else:
        if criteria == 0 and File.iat[i,11] < Low_resistance:
            Result.loc[Result.shape[0]] = File.iloc[i]
            criteria = 1
        else:
            continue

pd.DataFrame(Result).to_csv('Low_resistance'+'.txt', index=False)