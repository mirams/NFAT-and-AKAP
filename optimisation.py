import pints
from SpecialLogLikelihood import SpecialLogLikelihood
from model import NfatModel
import numpy as np
import matplotlib.pyplot as plt

# For reproducibility
np.random.seed(1)

log_like = SpecialLogLikelihood()

boundaries = pints.RectangularBoundaries(
    [1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9], [2, 2, 2, 2, 2, 2])

p_guess = np.array([1.65096678e-01, 3.31074825e-7, 1.30545310e-01, 2.12751554e-02,
                    1.99999983e+00, 2.03978544e-01])

print("Initial guess log likelihood = ", log_like(p_guess))

found_params, found_value = pints.optimise(
    log_like, p_guess, boundaries=boundaries, method=pints.CMAES)

print("Found params = ", found_params)
print("With score = ", log_like(found_params))

model = NfatModel()
times = np.linspace(-200, 360, 56001)
p = np.array([1, 1.2])
values, markers = model.simulate(
    np.concatenate([p, found_params, [0.0]]), times)

print("markers = ", markers)

# Plot the results
plt.figure()
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.plot(times, values)

a = values[:, 0]
ac = values[:, 1]
anc = values[:, 2]
ncp = values[:, 3]
nc = values[:, 4]
an = p[0] - (a + ac + anc)
nn = p[1] - (nc + an + anc + ncp)
percentage_akap_bound_to_nfat = 100*(an+anc)/p[0]
percentage_nfat_on_membrane = 100*(an+anc)/p[1]
plt.plot(times, nn)
plt.plot(times, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC",
           "NFATp (cytosol)", "NFAT (cytosol)", "NFAT (nucleus)", "AKAP:NFATp"])
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
