import pints
import pints.plot
from SpecialLogLikelihood import SpecialLogLikelihood
from model import NfatModel
import numpy as np
import matplotlib.pyplot as plt


import multiprocessing
multiprocessing.set_start_method('fork')

# For reproducibility
np.random.seed(0)

log_like = SpecialLogLikelihood()

boundaries = pints.RectangularBoundaries(
    [0, 0, 0, 0, 0], [2, 2, 2, 2, 2])

p_guess = np.array([2.31382813e-01, 3.14457334e-08,
                   2.64817597e-02, 2.10961798e-01, 1.38465592e-02])

print("Initial guess log likelihood = ", log_like(p_guess))

log_prior = pints.UniformLogPrior(
    [0, 0, 0, 0, 0],
    [2, 2, 2, 2, 2]
)

# Create a posterior log-likelihood (log(likelihood * prior))
log_posterior = pints.LogPosterior(log_like, log_prior)

nchains = 5
xs = pints.sample_initial_points(log_posterior, nchains, parallel=True)

# Create MCMC routine with four chains and run sampling.

mcmc = pints.MCMCController(log_posterior, nchains,
                            xs, method=pints.HaarioBardenetACMC)

# Add stopping criterion
mcmc.set_max_iterations(20000)

# Start adapting after 1000 iterations
mcmc.set_initial_phase_iterations(10000)

# Disable logging mode
mcmc.set_log_to_screen(True)

# Run!
print('Running...')
chains = mcmc.run()
print('Done!')

# Discard warm up
#chains = chains[:, 5000:, :]

# Look at distribution across all chains
pints.plot.pairwise(np.vstack(chains), kde=False)

# Show graphs
plt.show()

assert(False)
model = NfatModel()
times = np.linspace(-200, 60, 26001)
p = np.array([1, 1.2])
values, markers = model.simulate(np.concatenate([p, found_params]), times)

print("markers = ", markers)

# Plot the results
plt.figure()
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.plot(times, values)

a = values[:, 0]
ac = values[:, 1]
anc = values[:, 2]
nc = values[:, 3]
an = p[0] - (a + ac + anc)
nn = p[1] - (nc + an + anc)
percentage_akap_bound_to_nfat = 100*(an+anc)/p[0]
percentage_nfat_on_membrane = 100*(an+anc)/p[1]
plt.plot(times, an)
#plt.plot(times, nn)
plt.legend(["AKAP", "AKAP:CRAC", "AKAP:NFAT:CRAC",
           "NFAT (cytosol/nucleus)", "AKAP:NFAT"])
plt.show()

plt.figure()
plt.xlabel('Time')
plt.ylabel('Percentage (%)')
plt.plot(times, percentage_akap_bound_to_nfat)
plt.plot(times, percentage_nfat_on_membrane)
plt.plot(-1, 24, 'x')
plt.plot(30, 11, 'x')
plt.plot(-1, 26, 'x')
plt.plot(30, 17, 'x')
plt.legend(["AKAP bound to NFAT", "NFAT on membrane"])
plt.show()

akap_nfatp_times = np.genfromtxt(
    '../NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(0))
akap_nfatp_means = np.genfromtxt(
    '../NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(2))
akap_nfatp_sems = np.genfromtxt(
    '../NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(3))

plt.figure()
plt.xlabel('Time')
plt.ylabel('Proportion of AKAP:NFAT co-localised (% of value at t=0)')
plt.plot(akap_nfatp_times, 100*akap_nfatp_means, 'x')
plt.plot(times, 100*(an+anc)/np.max(an+anc))
plt.show()
