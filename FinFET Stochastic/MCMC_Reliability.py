import numpy as np
import pandas as pd
import sympy as sy
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from math import *
import sys

class MCMC:

    def MC_MX_sampler(data, burn_in, tol, num_cluster=None):  # Mixture Weibull sampler

        # Initialization
        num_data, dim_data = np.shape(data)
        cluster_index = [[] for i in range(num_cluster)]
        cluster_data = [[] for i in range(num_cluster)]
        if num_cluster == None:
            num_cluster = 2
        ini_list = np.random.multinomial(num_data, [1/num_cluster]*num_cluster, size=1)
        ini_index = np.arange(0, num_data, 1)
        for i in range(num_cluster):
            for j in range(ini_list[0, i]):
                t_ind =np.round(np.random.uniform(0, len(ini_index)-1)).astype(int)
                cluster_index[i].append(ini_index[t_ind])
                ini_index = np.delete(ini_index, t_ind)
            
        # Posterior Distribution sum(w_model*(theta*(s)**(alpha-1)*alpha*exp(-theta*(s)**alpha))
        theta_alpha = 0.001
        theta_beta = 0.001
        alpha_alpha = 1
        alpha_beta = 0.2
        w_fa = np.ones(num_cluster)
        w_model = np.random.dirichlet(w_fa)
        theta = np.random.gamma(theta_alpha, theta_beta, num_cluster)
        alpha = np.random.gamma(alpha_alpha, alpha_beta, num_cluster)

        # Update of P(z)
        for i in range(num_data):
            pz_list = np.array([w_model[k] * (theta[k]*alpha[k]*data[i]**(alpha[k]-1)) * exp(-theta[k]*data[i]**(alpha[k])) for k in range(num_cluster)])
            pz_list = np.cumsum(pz_list)/sum(pz_list)
            mv_pz, ind_pz = min([(mv_pz, ind_pz) for ind_pz, mv_pz in enumerate(abs(pz_list-np.random.uniform()))])
            cluster_data[ind_pz].append(data[i])

        # Update of w_model
        for i in range(num_cluster):
            w_model[i] = np.random.dirichlet(w_fa[i] + len(cluster_data[i]))

        # Update of theta
        for i in range(num_cluster):
            new_theta_alpha = len(cluster_data[i]) + theta_alpha
            new_theta_beta = theta_beta + sum([cluster_data[i][k]**(alpha[i]) for k in range(len(cluster_data[i]))])
            theta[i] = np.random.gamma(new_theta_alpha, new_theta_beta)

        # Update of alpha (sampled from Motroplis Hasting/Rejection support)
        for i in range(num_cluster):
            p_alpha = lambda x: x**(len(cluster_data[i]) + alpha_alpha - 1) * exp(x * sum([log(cluster_data[i][k]) for k in range(len(cluster_data[i]))]) - 
                theta[i] * sum([cluster_data[i][k]**(x) for k in range(len(cluster_data[i]))]) - x * alpha_beta)


        return



    def slicing_sampler(pdf, current_value, support=1e3): 
        
        P = pdf
        criteria = 1
        Uni_top = P(current_value)
        while criteria = 1:

            r = np.random.uniform(0, Uni_top)










        





def sampler1(data, samples=4, mu_init=.5, proposal_width=.5, plot=False, mu_prior_mu=0, mu_prior_sd=1.):
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
            plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_mu_ind, data, posterior_analytical)
        
        if accept:

            # Update position
            mu_current += mu_proposal
            mu_prososal = new_proposal
        
        posterior.append(mu_current)
        posterior_analytical.append(updates)
        
    return posterior


def sampler(data, samples=4, mu_init=.5, proposal_width=.5, plot=False, mu_prior_mu=0, mu_prior_sd=1.):
    mu_current = mu_init
    posterior = [mu_current]
    for i in range(samples):
        # suggest new position
        mu_proposal = norm(mu_current, proposal_width).rvs()

        # Compute likelihood by multiplying probabilities of each data point
        likelihood_current = norm(mu_current, 1).pdf(data).prod()
        likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
        likelihood_proposal = norm(num_cluster, 2).pdf(data).prod()
        
        # Compute prior probability of current and proposed mu
        # Define the pre-clustering factor as displayed by smapler 1        
        prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
        prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
        
        p_current = likelihood_current * prior_current
        p_proposal = likelihood_proposal * prior_proposal
        
        # Accept proposal?
        p_accept = p_proposal / p_current
        
        # Usually would include prior probability, which we neglect here for simplicity
        accept = np.random.rand() < p_accept
        reject -= new_proposal
        
        if plot:
            plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accept, posterior, i)
        
        if accept:
            # Update position on secondary powers 
            mu_current = mu_proposal
        
        posterior.append(mu_current)
        
    return posterior


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


# Function to display
def plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accepted, trace, i):
    from copy import copy
    trace = copy(trace)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, figsize=(16, 4))
    fig.suptitle('Iteration %i' % (i + 1))
    x = np.linspace(-3, 3, 5000)
    color = 'g' if accepted else 'r'
        
    # Plot prior
    prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
    prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
    prior = norm(mu_prior_mu, mu_prior_sd).pdf(x)
    ax1.plot(x, prior)
    ax1.plot([mu_current] * 2, [0, prior_current], marker='o', color='b')
    ax1.plot([mu_proposal] * 2, [0, prior_proposal], marker='o', color=color)
    ax1.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2), arrowprops=dict(arrowstyle="->", lw=2.))
    ax1.set(ylabel='Probability Density', title='current: prior(mu=%.2f) = %.2f\nproposal: prior(mu=%.2f) = %.2f' % (mu_current, prior_current, mu_proposal, prior_proposal))
    
    # Likelihood
    likelihood_current = norm(mu_current, 1).pdf(data).prod()
    likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
    y = norm(loc=mu_proposal, scale=1).pdf(x)
    sns.distplot(data, kde=False, norm_hist=True, ax=ax2)
    ax2.plot(x, y, color=color)
    ax2.axvline(mu_current, color='b', linestyle='--', label='mu_current')
    ax2.axvline(mu_proposal, color=color, linestyle='--', label='mu_proposal')
    #ax2.title('Proposal {}'.format('accepted' if accepted else 'rejected'))
    ax2.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    ax2.set(title='likelihood(mu=%.2f) = %.2f\nlikelihood(mu=%.2f) = %.2f' % (mu_current, 1e14*likelihood_current, mu_proposal, 1e14*likelihood_proposal))
    
    # Posterior
    posterior_analytical = calc_posterior_analytical(data, x, mu_prior_mu, mu_prior_sd)
    ax3.plot(x, posterior_analytical)
    posterior_current = calc_posterior_analytical(data, mu_current, mu_prior_mu, mu_prior_sd)
    posterior_proposal = calc_posterior_analytical(data, mu_proposal, mu_prior_mu, mu_prior_sd)
    ax3.plot([mu_current] * 2, [0, posterior_current], marker='o', color='b')
    ax3.plot([mu_proposal] * 2, [0, posterior_proposal], marker='o', color=color)
    ax3.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    #x3.set(title=r'prior x likelihood $\propto$ posterior')
    ax3.set(title='posterior(mu=%.2f) = %.5f\nposterior(mu=%.2f) = %.5f' % (mu_current, posterior_current, mu_proposal, posterior_proposal))
    
    if accepted:
        trace.append(mu_proposal)
    else:
        trace.append(mu_current)
    ax4.plot(trace)
    ax4.set(xlabel='iteration', ylabel='mu', title='trace')
    plt.tight_layout()
    #plt.legend()




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