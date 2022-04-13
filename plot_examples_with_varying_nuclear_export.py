from operator import delitem
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from model import NfatModel
from SpecialLogLikelihood import SpecialLogLikelihood
# Define the right-hand side of a system of ODEs


times = np.linspace(-200, 360, 56001)
# parameters
# p = [1, 1.2,
#      0.1195,
#      0.1003,
#      0.1089,
#      1.1357,
#      0.0453,
#      0]


# Parameters emphasising fit to time series data.
p = np.array([1, 1.2, 0.17743636, 0.0495394,  0.06270813,
             0.01552823, 0.34830096, 0.84645488, 0.1])

p = np.array([1, 1.2, 0.17743636, 0.0495394,  0.06270813,
             0.01552823, 0.34830096, 0.84645488, 0.5])

baseline_params = np.array(
    [1, 1.2, 0.17306501, 0.02632026, 0.07175625, 0.01612891, 0.37993729, 0.1])

# Parameters emphasising nuclear accumulation getting near 70%.
p = np.concatenate([baseline_params, [0]])  # NO NUCLEAR EXPORT


model = NfatModel()
values, markers = model.simulate(p, times)
print("markers = ", markers)


# Plot the results
a = values[:, 0]
ac = values[:, 1]
anc = values[:, 2]
ncp = values[:, 3]
nc = values[:, 4]
an = p[0] - (a + ac + anc)
nn = p[1] - (nc + an + anc + ncp)
prop_akap_bound_nfat_no_export = 100*(an+anc)/np.max(an+anc)

percentage_akap_bound_to_nfat = 100*(an+anc)/p[0]
percentage_nfat_on_membrane = 100*(an+anc)/p[1]
plt.plot(times, a)
plt.plot(times, ac)
plt.plot(times, anc)
plt.plot(times, ncp)
plt.plot(times, nc)
plt.plot(times, nn)
plt.plot(times, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC",
           "NFATp (cytosol)", "NFAT (cytosol)", "NFAT (nucleus)", "AKAP:NFATp"])
plt.show()

# Plot the results
plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Concentration (relative to total AKAP)')
plt.plot(times, a)
plt.plot(times, ac)
plt.plot(times, anc)
plt.plot(times, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC", "AKAP:NFATp"])
plt.plot([0, 0], [0, 1], '--', color='grey')
plt.plot([36, 36], [0, 1], '--', color='grey')
plt.xlim([-10, 200])
plt.show()


plt.figure()
plt.xlabel('Time (mins)')
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
    'data/NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(0))
akap_nfatp_means = np.genfromtxt(
    'data/NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(2))
akap_nfatp_means = akap_nfatp_means/np.mean(akap_nfatp_means[0:3])
# akap_nfatp_sems = np.genfromtxt(
#    '../NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(3))

plt.figure()
plt.xlabel('Time')
plt.ylabel('Proportion of AKAP:NFAT co-localised (% of mean value t<0)')
plt.plot(akap_nfatp_times, 100*akap_nfatp_means, 'x')
plt.plot(times, prop_akap_bound_nfat_no_export)
plt.xlim([-30, 150])
plt.plot([0, 0], [0, 110], '--', color='grey')
plt.plot([36, 36], [0, 110], '--', color='grey')
plt.show()


p = np.concatenate([baseline_params, [0.01]])  # SLOW NUCLEAR EXPORT

values, markers = model.simulate(p, times)
print("markers = ", markers)

# Plot the results
a = values[:, 0]
ac = values[:, 1]
anc = values[:, 2]
ncp = values[:, 3]
nc = values[:, 4]
an = p[0] - (a + ac + anc)
nn = p[1] - (nc + an + anc + ncp)
prop_akap_bound_nfat_slow_export = 100*(an+anc)/np.max(an+anc)

percentage_akap_bound_to_nfat = 100*(an+anc)/p[0]
percentage_nfat_on_membrane = 100*(an+anc)/p[1]
plt.plot(times, a)
plt.plot(times, ac)
plt.plot(times, anc)
plt.plot(times, ncp)
plt.plot(times, nc)
plt.plot(times, nn)
plt.plot(times, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC",
           "NFATp (cytosol)", "NFAT (cytosol)", "NFAT (nucleus)", "AKAP:NFATp"])
plt.show()

# Plot the results
plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Concentration (relative to total AKAP)')
plt.plot(times, a)
plt.plot(times, ac)
plt.plot(times, anc)
plt.plot(times, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC", "AKAP:NFATp"])
plt.plot([0, 0], [0, 1], '--', color='grey')
plt.plot([36, 36], [0, 1], '--', color='grey')
plt.xlim([-10, 200])
plt.show()


plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Percentage (%)')
plt.plot(times, percentage_akap_bound_to_nfat)
plt.plot(times, percentage_nfat_on_membrane)
plt.plot(-1, 24, 'x')
plt.plot(30, 11, 'x')
plt.plot(-1, 26, 'x')
plt.plot(30, 17, 'x')
plt.legend(["AKAP bound to NFAT", "NFAT on membrane"])
plt.show()


plt.figure()
plt.xlabel('Time')
plt.ylabel('Proportion of AKAP:NFAT co-localised (% of mean value t<0)')
plt.plot(akap_nfatp_times, 100*akap_nfatp_means, 'x')
plt.plot(times, prop_akap_bound_nfat_slow_export)
plt.xlim([-30, 150])
plt.plot([0, 0], [0, 110], '--', color='grey')
plt.plot([36, 36], [0, 110], '--', color='grey')
plt.show()

p = np.concatenate([baseline_params, [0.1]])   # IMPORT/EXPORT at same rate

values, markers = model.simulate(p, times)
print("markers = ", markers)

# Plot the results
a = values[:, 0]
ac = values[:, 1]
anc = values[:, 2]
ncp = values[:, 3]
nc = values[:, 4]
an = p[0] - (a + ac + anc)
nn = p[1] - (nc + an + anc + ncp)
prop_akap_bound_nfat_fast_export = 100*(an+anc)/np.max(an+anc)

percentage_akap_bound_to_nfat = 100*(an+anc)/p[0]
percentage_nfat_on_membrane = 100*(an+anc)/p[1]
plt.plot(times, a)
plt.plot(times, ac)
plt.plot(times, anc)
plt.plot(times, ncp)
plt.plot(times, nc)
plt.plot(times, nn)
plt.plot(times, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC",
           "NFATp (cytosol)", "NFAT (cytosol)", "NFAT (nucleus)", "AKAP:NFATp"])
plt.show()

# Plot the results
plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Concentration (relative to total AKAP)')
plt.plot(times, a)
plt.plot(times, ac)
plt.plot(times, anc)
plt.plot(times, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC", "AKAP:NFATp"])
plt.plot([0, 0], [0, 1], '--', color='grey')
plt.plot([36, 36], [0, 1], '--', color='grey')
plt.xlim([-10, 200])
plt.show()

plt.figure()
plt.xlabel('Time')
plt.ylabel('Proportion of AKAP:NFAT co-localised (% of mean value t<0)')
plt.plot(akap_nfatp_times, 100*akap_nfatp_means, 'x')
plt.plot(times, prop_akap_bound_nfat_no_export, color='orange')
plt.plot(times, prop_akap_bound_nfat_slow_export, color='red')
plt.plot(times, prop_akap_bound_nfat_fast_export, color='blue')
plt.legend(["Data", "Sim - no nuclear export",
           "Sim - slow nuclear export", "Sim - fast nuclear export"])
plt.xlim([-30, 150])
plt.plot([0, 0], [0, 110], '--', color='grey')
plt.plot([36, 36], [0, 110], '--', color='grey')
plt.show()
