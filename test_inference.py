import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import pints.toy

import multiprocessing
multiprocessing.set_start_method('fork')

# Load a forward model
model = pints.toy.LogisticModel()

# Create some toy data
real_parameters = [0.015, 500]
times = np.linspace(0, 1000, 1000)
org_values = model.simulate(real_parameters, times)

# Add noise
noise = 10
values = org_values + np.random.normal(0, noise, org_values.shape)
real_parameters = np.array(real_parameters + [noise])

# Get properties of the noise sample
noise_sample_mean = np.mean(values - org_values)
noise_sample_std = np.std(values - org_values)

# Create an object with links to the model and time series
problem = pints.SingleOutputProblem(model, times, values)

# Create a log-likelihood function (adds an extra parameter!)
log_likelihood = pints.GaussianLogLikelihood(problem)

# Create a uniform prior over both the parameters and the new noise variable
log_prior = pints.UniformLogPrior(
    [0.01, 400, noise*0.1],
    [0.02, 600, noise*100]
)

# Create a posterior log-likelihood (log(likelihood * prior))
log_posterior = pints.LogPosterior(log_likelihood, log_prior)

# Select random starting location for each chain: this method will sample from the log_prior distribution since a pints.LogPosterior has been supplied as the objective. (Note, this can be overridden by supplying a random sampling function as an input.)

nchains = 4
xs = pints.sample_initial_points(log_posterior, nchains, parallel=True)

# Create MCMC routine with four chains and run sampling.

mcmc = pints.MCMCController(log_posterior, nchains,
                            xs, method=pints.HaarioBardenetACMC)

# Add stopping criterion
mcmc.set_max_iterations(10000)

# Start adapting after 1000 iterations
mcmc.set_initial_phase_iterations(1000)

# Disable logging mode
mcmc.set_log_to_screen(False)

# Run!
print('Running...')
chains = mcmc.run()
print('Done!')

# Discard warm up
chains = chains[:, 5000:, :]

# Look at distribution across all chains
pints.plot.pairwise(np.vstack(chains), kde=False)

# Show graphs
plt.show()

results = pints.MCMCSummary(
    chains=chains, time=mcmc.time(), parameter_names=['r', 'k', 'sigma'])
print(results)
