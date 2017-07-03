import numpy as np
import pandas as pd
import sympy as sy
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import statsmodels.tsa.api as sm
from math import *

class MCMC:

    def MCMC_MX_sampler(data, burn_in, test, tol, num_cluster=None, thinning_gap=100):  # Mixture Weibull sampler

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
        theta_beta = 0.001
        alpha_alpha = 1
        alpha_beta = 0.2
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

                pz_list = np.array([w_model[k] * (theta[k]*alpha[k]*data[i]**(alpha[k]-1)) * exp(-theta[k]*data[i]**(alpha[k])) for k in range(num_cluster)])

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
                new_theta[i] = np.random.gamma(new_theta_alpha, new_theta_beta)

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




    def MCMC_MC_sampler(data, burn_in, test, tol, num_cluster=None):  # Mixture Weibull sampler

        # Initialization
        num_data = len(data)
        if num_cluster == None:
            num_cluster = 2
        count= 0
            
        # Posterior Distribution sum(w_model*(theta*(s)**(alpha-1)*alpha*exp(-theta*(s)**alpha))
        theta_alpha = 0.01
        theta_beta = 0.001
        alpha_alpha = 1
        alpha_beta = 0.2
        w_fa = np.ones(num_cluster)
        origin_w_fa = w_fa[:]
        w_model = np.random.dirichlet(w_fa)
        theta = np.random.gamma(theta_alpha, theta_beta, num_cluster)
        alpha = np.random.gamma(alpha_alpha, alpha_beta, num_cluster)
        likelihood_record = []
        w_record = []
        theta_record = []
        alpha_record = []

        while count <= (burn_in + test):

            count += 1
            log_likelihood = 0
            cluster_data = [[] for i in range(num_cluster)]
            new_theta = [[] for i in range(num_cluster)]
            new_alpha = [[] for i in range(num_cluster)]

            # Update of P(z)
            #import pdb; pdb.set_trace()
            for i in range(num_data):

                pz_list = np.array([w_model[k] * (theta[k]*alpha[k]*data[i]**(alpha[k]-1)) * exp(-theta[k]*data[i]**(alpha[k])) for k in range(num_cluster)])

                mv_pz, ind_pz = max([(mv_pz, ind_pz) for ind_pz, mv_pz in enumerate(pz_list)])
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
                new_theta[i] = np.random.gamma(new_theta_alpha, new_theta_beta)

            theta = new_theta[:]

            # Update of alpha (sampled from Motroplis Hasting/Rejection support)
            #import pdb; pdb.set_trace()
            for i in range(num_cluster):

                p_alpha = lambda x: x**(len(cluster_data[i]) + alpha_alpha - 1) * exp(x * sum([log(cluster_data[i][k]) for k in range(len(cluster_data[i]))]) - 
                    theta[i] * sum([cluster_data[i][k]**(x) for k in range(len(cluster_data[i]))]) - x * alpha_beta)

                new_alpha[i] = MCMC.slice_sampler(p_alpha, alpha[i], left_bound=0)
            
            alpha = new_alpha[:]


            # Record the sampling after burn-in
            if count >= burn_in:

                w_record.append(w_model[:])
                theta_record.append(theta[:])
                alpha_record.append(alpha[:])

            likelihood_record.append(log_likelihood)
            

        return w_record, theta_record, alpha_record, likelihood_record





    def slice_sampler(pdf, current_value, step=10, left_bound=None, right_bound=None): 

        P = pdf
        criteria = 1
        Uni_top = P(current_value)
        count = 0
        left = -step
        right = step
        ad_l_bound = left_bound
        ad_r_bound = right_bound
        while criteria == 1:

            count += 1
            # find left boundary
            left_out = MCMC.find_boundary(P, current_value, Uni_top, left, l_bound=ad_l_bound, r_bound=ad_r_bound)
            # find right boundary
            right_out = MCMC.find_boundary(P, current_value, Uni_top, right, l_bound=ad_l_bound, r_bound=ad_r_bound)

            n_point = np.random.uniform(0,1) * (right_out - left_out) + left_out

            if P(n_point) >= Uni_top:
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
                    direction *= 2

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
                    direction *= 2

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
                    direction *= 2

        else:

            while criteria == 1:

                n_point = s_point + direction
                if pdf(n_point) <= support:
                    out = n_point
                    criteria = 0
                    break
                else:
                    direction *= 2

        return out




    def data_preprocessing(data, shift=1e-9):

        data = data.astype(np.float32, copy = False)
        data = data[~np.isnan(data)]
        min_value = min(data)
        scale = max(data) - min(data)

        data = (data -min_value) /  scale + shift
        
        return data, scale, min_value




    def model_reconstruction(w_record, theta_record, alpha_record, raw_data, scale, min_value, shift=1e-9):

        num_cluster = np.size(w_record, 0)

        data = np.sort(raw_data, axis = 0)
        data_length = len(data)
        raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
        plot_p = np.log(-np.log(1 - raw_probability))

        w = np.mean(w_record, axis=1)
        theta = np.mean(theta_record, axis=1)
        alpha = np.mean(alpha_record, axis=1)

        w_std = np.std(w_record, axis=1)
        theta_std = np.std(theta_record, axis=1)
        alpha_std = np.std(alpha_record, axis=1)

        summary = [(w, w_std), (theta, theta_std), (alpha, alpha_std)]

        fiiting_range = np.logspace(np.log(shift), np.log(1+shift), 200, base=e)
        real_fit_range = (fiiting_range - shift) * scale + min_value

        fitting_probability = np.array([1 - sum([ w[k] * exp(- theta[k] * t ** (alpha[k])) for k in range(num_cluster)]) for t in fiiting_range])
        plot_fit_p = np.log(-np.log(1 - fitting_probability))



        return real_fit_range, plot_fit_p, data, plot_p, summary


    def model_reconstruction_1(w_record, theta_record, alpha_record, raw_data):

        num_cluster = np.size(w_record, 0)

        data = np.sort(raw_data, axis = 0)
        data_length = len(data)
        raw_probability = np.array([(i - 0.3) / data_length for i in range(1, data_length+1)])
        plot_p = np.log(-np.log(1 - raw_probability))

        w = np.mean(w_record, axis=1)
        theta = np.mean(theta_record, axis=1)
        alpha = np.mean(alpha_record, axis=1)

        w_std = np.std(w_record, axis=1)
        theta_std = np.std(theta_record, axis=1)
        alpha_std = np.std(alpha_record, axis=1)

        summary = [(w, w_std), (theta, theta_std), (alpha, alpha_std)]

        fiiting_range = np.logspace(np.log(min(data)), np.log(max(data)), 200, base=e)

        fitting_probability = np.array([1 - sum([ w[k] * exp(- theta[k] * t ** (alpha[k])) for k in range(num_cluster)]) for t in fiiting_range])
        plot_fit_p = np.log(-np.log(1 - fitting_probability))



        return fiiting_range, plot_fit_p, data, plot_p, summary








def BFGS(expr, initial, tol, alpha = 1, mix_iter = 50):

    var_list = list(initial.keys())
    temp_value = list(initial.values())
    object_func = sy.lambdify((var_list,), expr, 'numpy', dummify=False)    
    gradient = [sy.diff(expr, i) for i in var_list]
    gradient_func = sy.lambdify((var_list,), gradient, 'numpy', dummify=False)

    error = 1
    I_ini = np.identity(len(var_list))
    hessian_inv = I_ini
    count = 0

    #import pdb; pdb.set_trace()


    while error >= tol:

        g = np.array(gradient_func(temp_value))
        p = -1.0 * np.dot(hessian_inv, g)
        p = p / np.dot(p,p)
        step = step_checking(object_func, gradient_func, p, temp_value, 1e3, tol)
        s = step * p
        if step == 0:
            count += 1
            break

        new_value = temp_value + s
        new_g = np.array(gradient_func(new_value))
        y = new_g - g
        # Define secondary qusai newtown updates
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



def step_checking(f, g, d, x, ini_alpha, tol, mixiter=50):
    rho = 0.2
    sigma = 1/6
    step_criteria = 0
    alpha = ini_alpha
    final_a = 1
    count = 0

    #import pdb; pdb.set_trace()

    while step_criteria == 0:
        count += 1
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
            alpha = alpha*rho/sigma

        if count > mixiter:
            break


    return final_a





Data_file = pd.ExcelFile(r'C:\Users\Mason\Documents\Project\Matlab Project\Clustering data processing\FinFET\125 TDDB FR0709_position.xlsx')  # revise path
#Data_file = pd.ExcelFile(r'C:\Users\Sen\Desktop\1.xlsx')
p_data = Data_file.parse('Sheet1', index_row = None, header = None)
p_data.drop(p_data.columns[[0]], axis = 0, inplace  =True)  # drop first row
p_data = p_data.iloc[:,:].values
data = p_data[:,0]    # data from which column
data = data.astype(np.float32, copy = False)
data = data[~np.isnan(data)]
Data, scale, min_value = MCMC.data_preprocessing(data)


#np.seterr(divide='ignore', invalid='ignore', over='ignore')

set_burn_in=1e5
set_test=1e5
num_cluster_set=3
w_record, theta_record, alpha_record, likelihood_record, Autocorrelation_w, Autocorrelation_theta, Autocorrelation_alpha = MCMC.MCMC_MX_sampler(Data, burn_in=set_burn_in, test=set_test, tol=1e-9, num_cluster=num_cluster_set, thinning_gap=1000)
real_fit_range, plot_fit_p, real_data, plot_p, summary = MCMC.model_reconstruction(w_record, theta_record, alpha_record, data, scale, min_value)
#real_fit_range, plot_fit_p, real_data, plot_p, summary = MCMC.model_reconstruction_1(w_record, theta_record, alpha_record, data)

#import pdb; pdb.set_trace()
plt.interactive(True)
f1, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
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
ax2.plot(real_data, plot_p, 'bo', markersize=10) 
ax2.plot(real_fit_range, plot_fit_p, 'b-', linewidth=2)

print(summary)



# plt.interactive(True)
# plt.plot(range(int(set_burn_in)+101), likelihood_record, 'bo', markersize=10) 
# plt.xlabel(r'Iteration',{'fontname':'Times New Roman','fontsize':18})
# plt.ylabel(r'Loglikelihood',{'fontname':'Times New Roman','fontsize':18})

