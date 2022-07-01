import numpy as np
import matplotlib.pyplot as plt
from model import NfatModel

times = np.linspace(-200, 360, 56001)

# Old params from original optimsation (see change on lines 42-44 of specialloglikelihood)
# baseline_params = np.array(
#    [1, 1.2, 0.17299323, 0.0273452,  0.07162728, 0.01610738, 0.38311162, 0.1])

baseline_params = np.array(
    [1, 2, 0.14611337, 0.02547575, 0.07700975, 0.01420659, 0.59655061, 0.1])

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
plt.plot(times+10, a)
plt.plot(times+10, ac)
plt.plot(times+10, anc)
plt.plot(times+10, ncp)
plt.plot(times+10, nc)
plt.plot(times+10, nn)
plt.plot(times+10, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC",
           "NFATp (cytosol)", "NFAT (cytosol)", "NFAT (nucleus)", "AKAP:NFATp"])
plt.show()

# Plot the results
plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Concentration (relative to total AKAP)')
plt.plot(times+10, a)
plt.plot(times+10, ac)
plt.plot(times+10, anc)
plt.plot(times+10, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC", "AKAP:NFATp"])
plt.plot([10, 10], [0, 1], '--', color='grey')
plt.plot([46, 46], [0, 1], '--', color='grey')
plt.xlim([0, 210])
plt.show()


plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Percentage (%)')
plt.plot(times+10, percentage_akap_bound_to_nfat)
plt.plot(times+10, percentage_nfat_on_membrane)
plt.plot(9, 12, 'x')
plt.plot(40, 5.6, 'x')
plt.plot(9, 26, 'x')
plt.plot(40, 17, 'x')
plt.legend(["AKAP bound to NFAT", "NFAT on membrane"])
plt.show()


akap_nfatp_times = np.genfromtxt(
    'data/NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(0))
akap_nfatp_means = np.genfromtxt(
    'data/NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(2))
akap_nfatp_means = akap_nfatp_means/np.mean(akap_nfatp_means[0:3])
akap_nfatp_sems = np.genfromtxt(
    'data/NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(3))/np.mean(akap_nfatp_means[0:3])

plt.figure()
plt.xlabel('Time')
plt.ylabel('Proportion of AKAP:NFAT co-localised (% of mean value t<0)')
plt.errorbar(akap_nfatp_times+10, 100*akap_nfatp_means,
             yerr=100*akap_nfatp_sems, capsize=1, lw=1, ls='', marker='x')
plt.plot(times+10, prop_akap_bound_nfat_no_export)
plt.xlim([0, 160])
plt.plot([10, 10], [0, 110], '--', color='grey')
plt.plot([46, 46], [0, 110], '--', color='grey')
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
plt.plot(times+10, a)
plt.plot(times+10, ac)
plt.plot(times+10, anc)
plt.plot(times+10, ncp)
plt.plot(times+10, nc)
plt.plot(times+10, nn)
plt.plot(times+10, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC",
           "NFATp (cytosol)", "NFAT (cytosol)", "NFAT (nucleus)", "AKAP:NFATp"])
plt.show()

# Plot the results
plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Concentration (relative to total AKAP)')
plt.plot(times+10, a)
plt.plot(times+10, ac)
plt.plot(times+10, anc)
plt.plot(times+10, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC", "AKAP:NFATp"])
plt.xlim([0, 160])
plt.plot([10, 10], [0, 1], '--', color='grey')
plt.plot([46, 46], [0, 1], '--', color='grey')
plt.show()


plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Percentage (%)')
plt.plot(times+10, percentage_akap_bound_to_nfat)
plt.plot(times+10, percentage_nfat_on_membrane)
plt.plot(-1, 12, 'x')
plt.plot(30, 5.6, 'x')
plt.plot(-1, 26, 'x')
plt.plot(30, 17, 'x')
plt.legend(["AKAP bound to NFAT", "NFAT on membrane"])
plt.show()


plt.figure()
plt.xlabel('Time')
plt.ylabel('Proportion of AKAP:NFAT co-localised (% of mean value t<0)')
plt.errorbar(akap_nfatp_times+10, 100*akap_nfatp_means,
             yerr=100*akap_nfatp_sems, capsize=1, lw=1, ls='', marker='x')
plt.plot(times+10, prop_akap_bound_nfat_slow_export)
plt.xlim([-20, 160])
plt.plot([10, 10], [0, 110], '--', color='grey')
plt.plot([46, 46], [0, 110], '--', color='grey')
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
plt.plot(times+10, a)
plt.plot(times+10, ac)
plt.plot(times+10, anc)
plt.plot(times+10, ncp)
plt.plot(times+10, nc)
plt.plot(times+10, nn)
plt.plot(times+10, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC",
           "NFATp (cytosol)", "NFAT (cytosol)", "NFAT (nucleus)", "AKAP:NFATp"])
plt.show()

# Plot the results
plt.figure()
plt.xlabel('Time (mins)')
plt.ylabel('Concentration (relative to total AKAP)')
plt.plot(times+10, a)
plt.plot(times+10, ac)
plt.plot(times+10, anc)
plt.plot(times+10, an)
plt.legend(["AKAP (free)", "AKAP:CRAC", "AKAP:NFATp:CRAC", "AKAP:NFATp"])
plt.xlim([0, 160])
plt.plot([10, 10], [0, 1], '--', color='grey')
plt.plot([46, 46], [0, 1], '--', color='grey')
plt.show()

plt.figure()
plt.xlabel('Time (minutes since Thap addition)')
plt.ylabel('Proportion of AKAP co-localised with NFAT (% of mean value t<0)')
plt.errorbar(akap_nfatp_times+10, 100*akap_nfatp_means,
             yerr=100*akap_nfatp_sems, capsize=1, lw=1, ls='', marker='x')
plt.plot(times+10, prop_akap_bound_nfat_no_export, color='orange')
plt.plot(times+10, prop_akap_bound_nfat_slow_export, color='red')
plt.plot(times+10, prop_akap_bound_nfat_fast_export, color='blue')
plt.legend(["Sim - no nuclear export",
           "Sim - slow nuclear export", "Sim - fast nuclear export", "Data +/- s.e.m."])
plt.xlim([-20, 160])
plt.plot([10, 10], [0, 110], '--', color='grey')
plt.plot([46, 46], [0, 110], '--', color='grey')
plt.show()
