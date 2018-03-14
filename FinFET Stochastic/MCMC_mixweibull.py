import numpy as np
import pandas as pd
import sympy as sy
from multiprocessing import Pool
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import statsmodels.tsa.api as sm
import pymc3 as pm
from math import *

class MCMC:

    def MCMC_MW_sampler(data, burn_in, test, tol, num_cluster=None, thinning_gap=10):  # Mixture Weibull sampler

        # Initialization
        num_data = len(data)
        if num_cluster == None:
            num_cluster = 2
        count = 0
        thin_count = 0
        acm_length = int(np.round(burn_in + test)/thinning_gap)
            
        # Posterior Distribution sum(w_model*(theta*(s)**(alpha-1)*alpha*exp(-theta*(s)**alpha))
        # Pramaters
        theta_alpha = 0.01
        theta_beta = 1 / 0.01
        alpha_alpha = 1
        alpha_beta = 1 / 0.2
        w_fa = np.ones(num_cluster)
        origin_w_fa = w_fa[:]
        w_model = np.random.dirichlet(w_fa)
        theta = np.random.gamma(theta_alpha, theta_beta, num_cluster)
        alpha = np.random.gamma(alpha_alpha, alpha_beta, num_cluster)
        # Samples
        likelihood_record = np.zeros(int(burn_in + test))
        w_record = np.zeros((num_cluster,int(test)))
        theta_record = np.zeros((num_cluster,int(test)))
        alpha_record = np.zeros((num_cluster,int(test)))
        # Autocorrelation Check
        w_acm = np.zeros((num_cluster, acm_length))
        theta_acm = np.zeros((num_cluster, acm_length))
        alpha_acm = np.zeros((num_cluster, acm_length))

        while count <= (burn_in + test - 1):

            count += 1
            log_likelihood = 0
            cluster_data = [[] for i in range(num_cluster)]
            new_theta = [[] for i in range(num_cluster)]
            new_alpha = [[] for i in range(num_cluster)]

            # Update of P(z)
            #import pdb; pdb.set_trace()
            for i in range(num_data):

                pz_list = np.array([w_model[k] * (theta[k] * alpha[k] * (data[i]**(alpha[k]-1))) * exp( - theta[k] * (data[i]**(alpha[k]))) for k in range(num_cluster)])

                total_pz_list = np.cumsum(pz_list)
                total_pz_list = total_pz_list - np.random.uniform(0,1) * total_pz_list[-1]
                total_pz_list[total_pz_list < 0] = inf
                mv_pz, ind_pz = min([(mv_pz, ind_pz) for ind_pz, mv_pz in enumerate(total_pz_list)])
                cluster_data[ind_pz].append(data[i])

                if pz_list[ind_pz] == 0:
                    log_likelihood += -1/tol
                else:
                    log_likelihood += log(pz_list[ind_pz])

            # Update of w_model
            #import pdb; pdb.set_trace()
            for i in range(num_cluster):

                w_fa[i] = origin_w_fa[i] + len(cluster_data[i])

            w_model = np.random.dirichlet(w_fa)


            # Update of theta
            #import pdb; pdb.set_trace()
            for i in range(num_cluster):

                new_theta_alpha = len(cluster_data[i]) + theta_alpha
                new_theta_beta = theta_beta + sum([cluster_data[i][k]**(alpha[i]) for k in range(len(cluster_data[i]))])
                new_theta[i] = np.random.gamma(new_theta_alpha, 1 / new_theta_beta)

            theta = new_theta[:]

            # Update of alpha (sampled from Motroplis Hasting/Rejection support)
            #import pdb; pdb.set_trace()
            for i in range(num_cluster):

                p_alpha = lambda x: x**(len(cluster_data[i]) + alpha_alpha - 1) * exp(x * sum([log(cluster_data[i][k]) for k in range(len(cluster_data[i]))]) - 
                    theta[i] * sum([cluster_data[i][k]**(x) for k in range(len(cluster_data[i]))]) - x * alpha_beta)

                new_alpha[i] = MCMC.slice_sampler(p_alpha, alpha[i], left_bound=0)
            
            alpha = new_alpha[:]


            # Record the sampling after burn-in
            if count > burn_in:
                w_record[:,int(count - burn_in - 1)] = w_model[:]
                theta_record[:,int(count - burn_in - 1)] = theta[:]
                alpha_record[:,int(count - burn_in - 1)] = alpha[:]

            likelihood_record[count - 1] = log_likelihood

            # Thinning the Markov Chain and check the convergency
            if count == (thin_count + 1) * thinning_gap:
                w_acm[:,thin_count] = w_model[:]
                theta_acm[:,thin_count] = theta[:]
                alpha_acm[:,thin_count] = alpha[:]
                thin_count += 1

        Autocorrelation_w = [sm.stattools.acf(w_acm[i], unbiased=True, nlags=int(np.size(w_acm)/2), fft=True) for i in range(num_cluster)]
        Autocorrelation_theta = [sm.stattools.acf(theta_acm[i], unbiased=True, nlags=int(np.size(w_acm)/2), fft=True) for i in range(num_cluster)]
        Autocorrelation_alpha = [sm.stattools.acf(alpha_acm[i], unbiased=True, nlags=int(np.size(w_acm)/2), fft=True) for i in range(num_cluster)]

      
        return w_record, theta_record, alpha_record, likelihood_record, Autocorrelation_w, Autocorrelation_theta, Autocorrelation_alpha


    def slice_sampler(pdf, current_value, step=0.2, left_bound=None, right_bound=None): 

        P = pdf
        criteria = 1
        Uni_top = P(current_value)
        y = Uni_top * np.random.uniform(0,1)
        count = 0
        left = -step
        right = step
        ad_l_bound = left_bound
        ad_r_bound = right_bound
        while criteria == 1:

            count += 1
            # find left boundary
            left_out = MCMC.find_boundary(P, current_value, y, left, l_bound=ad_l_bound, r_bound=ad_r_bound)
            # find right boundary
            right_out = MCMC.find_boundary(P, current_value, y, right, l_bound=ad_l_bound, r_bound=ad_r_bound)

            n_point = np.random.uniform(0,1) * (right_out - left_out) + left_out

            if P(n_point) >= y:
                new_sample = n_point
                criteria = 0
                break
            elif n_point >= current_value:

                ad_r_bound = n_point
            else:

                ad_l_bound = n_point

        return new_sample


    def find_boundary(pdf, s_point, support, direction, l_bound=None, r_bound=None):

        criteria = 1

        fix_point = s_point
        count = 0

        if l_bound == None and r_bound != None:

            while criteria == 1:

                n_point = s_point + direction
                if n_point >= r_bound:
                    out = r_bound
                    criteria = 0
                    break
                elif pdf(n_point) <= support:
                    out = n_point
                    criteria = 0
                    break
                else:
                    count += 1
                    s_point = n_point
                    if count > 2000:
                       criteria = 0 
                       out = fix_point

        elif l_bound != None and r_bound == None:

            while criteria == 1:

                n_point = s_point + direction
                if n_point <= l_bound:
                    out = l_bound
                    criteria = 0
                    break
                elif pdf(n_point) <= support:
                    out = n_point
                    criteria = 0
                    break
                else:
                    count += 1
                    s_point = n_point
                    if count > 2000:
                       criteria = 0
                       out = fix_point

        elif l_bound != None and r_bound != None:

            while criteria == 1:

                n_point = s_point + direction
                if n_point <= l_bound:
                    out = l_bound
                    criteria = 0
                    break
                elif n_point >= r_bound:
                    out = r_bound
                    criteria = 0
                    break                    
                elif pdf(n_point) <= support:
                    out = n_point
                    criteria = 0
                    break
                else:
                    count += 1
                    s_point = n_point
                    if count > 2000:
                       criteria = 0
                       out = fix_point

        else:

            while criteria == 1:

                n_point = s_point + direction
                if pdf(n_point) <= support:
                    out = n_point
                    criteria = 0
                    break
                else:
                    count += 1
                    s_point = n_point
                    if count > 2000:
                       criteria = 0
                       out = fix_point

        return out


    def data_preprocessing(data, shift=1e-9):

        data = data.astype(np.float32, copy = False)
        data = data[~np.isnan(data)]
        min_value = min(data)
        scale = max(data) 

        data = (data) / scale 
        
        return data, scale, min_value


    def model_reconstruction(w_record, theta_record, alpha_record, raw_data):

        num_cluster = np.size(w_record, 0)

        data = np.sort(raw_data, axis = 0)
        data_length = len(data)
        raw_probability = np.array([(i - 0.3) / (data_length + 0.4) for i in range(1, data_length+1)])
        plot_p = np.log(-np.log(1 - raw_probability))

        w = np.mean(w_record, axis=1)
        theta = np.mean(theta_record, axis=1)
        alpha = np.mean(alpha_record, axis=1)

        
        w_std = np.std(w_record, axis=1)
        theta_std = np.std(theta_record, axis=1)
        alpha_std = np.std(alpha_record, axis=1)

        summary = [(w, w_std), (theta, theta_std), (alpha, alpha_std)]

        fiiting_range = np.logspace(np.log(min(raw_data)), np.log(max(raw_data)), 200, base=e)

        fitting_probability = np.array([1 - sum([ w[k] * exp(- theta[k] * t ** (alpha[k])) for k in range(num_cluster)]) for t in fiiting_range])
        plot_fit_p = np.log(-np.log(1 - fitting_probability))


        return fiiting_range, plot_fit_p, data, plot_p, summary


Data_file = pd.ExcelFile(r'C:\Users\Mason\Documents\Project\Matlab Project\Clustering data processing\FinFET\MCMC.xlsx')  # revise path
#Data_file = pd.ExcelFile(r'C:\Users\Sen\Desktop\1.xlsx') 
p_data = Data_file.parse('Sheet1', index_row = None, header = None)
p_data.drop(p_data.columns[[0]], axis = 0, inplace  =True)  # drop first row
p_data = p_data.iloc[:,:].values
data = p_data[:,3]    # data from which column
data = data.astype(np.float32, copy = False)
data = data[~np.isnan(data)]
Data, scale, min_value = MCMC.data_preprocessing(data)

#MCMC Setting

set_burn_in = 2e4
set_test = 1e4
num_cluster_set = 4
set_thinning_gap = 100
set_tol = 1e-9


#Multiprocess and BIC Checking

number_BIC = 1
number_cores = 3


input_list = zip([data for i in range(int(number_BIC))], [set_burn_in for i in range(int(number_BIC))], [set_test for i in range(int(number_BIC))], 
    [set_tol for i in range(int(number_BIC))], [num_cluster_set for i in range(int(number_BIC))], [set_thinning_gap for i in range(int(number_BIC))])

#import pdb; pdb.set_trace()

if __name__ == "__main__":
    #w_record, theta_record, alpha_record, likelihood_record, Autocorrelation_w, Autocorrelation_theta, Autocorrelation_alpha = MCMC.MCMC_MW_sampler(data, burn_in=set_burn_in, test=set_test, tol=1e-9, num_cluster=num_cluster_set, thinning_gap=100)
    with Pool(number_cores) as pool:
        Results = pool.starmap(MCMC.MCMC_MW_sampler, input_list)
        pool.close()
        pool.join()

#BIC
BIC = [np.log(len(data)) * num_cluster_set * 3 - 2 * np.mean(Results[i][4]) for i in range(int(number_BIC))]
mv_BIC, ind_BIC = min([(mv_BIC, ind_BIC) for ind_BIC, mv_BIC in enumerate(BIC)])

w_record, theta_record, alpha_record, likelihood_record, Autocorrelation_w, Autocorrelation_theta, Autocorrelation_alpha = Results[ind_BIC]

fiiting_range, plot_fit_p, real_data, plot_p, summary = MCMC.model_reconstruction(w_record, theta_record, alpha_record, data)

#Plot and Visulization

plt.interactive(True)
f1, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
for i in range(num_cluster_set):
    ax1.plot(range(int(np.size(Autocorrelation_w[i]))), Autocorrelation_w[i],'b-')
    ax2.plot(range(int(np.size(Autocorrelation_theta[i]))), Autocorrelation_theta[i],'r-')
    ax3.plot(range(int(np.size(Autocorrelation_alpha[i]))), Autocorrelation_alpha[i],'g-')
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
ax1.set_title('Convergence of Decision Tree')
f1.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f1.axes[:-1]], visible=False)

f2 = plt.figure()
ax2 = f2.add_subplot(111)
ax2.plot(real_data, plot_p, 'ko', markersize=6) 
ax2.plot(fiiting_range, plot_fit_p, 'g-', linewidth=4)


ax2.set_xscale('log')   
ax2.set_xlabel(r'Time to Failure (s)',{'fontname':'Times New Roman','fontsize':18})
ax2.set_ylabel(r'ln(-ln(1-F))',{'fontname':'Times New Roman','fontsize':18})

f3 = plt.figure()
ax3 = f3.add_subplot(111)
ax3.plot(range(len(likelihood_record)), -  likelihood_record, 'b-', markersize=10)

ax3.set_xscale('log')


print(summary)
print(np.mean(likelihood_record[int(set_burn_in):int(set_burn_in + set_test -1)], axis=0))


