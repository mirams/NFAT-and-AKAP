import numpy as np
import pints
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the right-hand side of a system of ODEs


class NfatModel(pints.ForwardModel):

    def simulate(self, parameters, times):

        def rhs(y, t, params):

            assert(len(params) == 9)

            a_total = params[0]  # Total AKAP
            n_total = params[1]  # Total NFAT
            k1 = params[2]  # ANp association with CRAC
            k2 = params[3]  # ANp dissociation from ANC
            k3 = params[4]  # NFATp dissociation from AKAP
            k4 = params[5]  # Ncp to AKAP binding (whether to A or AC)
            k5 = params[6]  # ANC -> Nc when CRAC ON
            k6 = params[7]  # Nc nuclear import
            k7 = params[8]  # Ncp nuclear export

            if t <= 0 or t > 36:
                k1 = 0
                k5 = 0

            a = y[0]
            ac = y[1]
            anc = y[2]
            ncp = y[3]
            nc = y[4]
            an = a_total - (a + ac + anc)
            nn = n_total - (ncp + an + anc + nc)

            dy = np.zeros(5)

            dy[0] = - k1*a + k2*ac - k4*a*ncp + k3*an
            dy[1] = + k1*a - k2*ac + (k3+k5)*anc - k4*ac*ncp
            dy[2] = + k1*an - k2*anc + k4*ac*ncp - (k3+k5)*anc
            dy[3] = + k3*anc - k4*(ac+a)*ncp + k7*nn + k3*an
            dy[4] = + k5*anc - k6*nc
            return dy

        # Run a simulation
        y0 = [parameters[0], 0, 0, parameters[1], 0]      # initial conditions

        # Call odeint, with the parameters wrapped in a tuple
        values = odeint(rhs, y0, times, (parameters,), printmessg=1, hmax=10)

        # The markers we want from the outputs:
        markers = np.zeros(6)
        a = values[:, 0]
        ac = values[:, 1]
        anc = values[:, 2]
        # ncp = values[: , 3] # unused
        # nc = values[:,4] # Unused
        an = parameters[0] - (a + ac + anc)
        # nn = parameters[1] - (nc + an + anc + ncp) # unused
        percentage_akap_bound_to_nfat = 100*(an+anc)/parameters[0]
        percentage_nfat_on_membrane = 100*(an+anc)/parameters[1]

        idx_time_zero = np.argmax(times >= 0)
        idx_steady_on = np.argmax(times >= 36)

        markers[0] = percentage_nfat_on_membrane[idx_time_zero]
        markers[1] = percentage_nfat_on_membrane[idx_steady_on]
        markers[2] = percentage_akap_bound_to_nfat[idx_time_zero]
        markers[3] = percentage_akap_bound_to_nfat[idx_steady_on]

        # Rough and ready half-life for drop
        idx3 = np.where(
            percentage_nfat_on_membrane[idx_time_zero:idx_steady_on] < markers[1] + 0.5*(markers[0] - markers[1]))[0]
        if idx3.size != 0:
            markers[4] = times[idx_time_zero+idx3[0]]
        else:
            markers[4] = 1e4

        # and recovery
        idx4 = np.where(
            percentage_nfat_on_membrane[idx_steady_on:] > markers[1] + 0.5*(markers[0] - markers[1]))[0]
        if idx4.size != 0:
            markers[5] = times[idx_steady_on+idx4[0]]-times[idx_steady_on]
        else:
            markers[5] = 1e4

        return values, markers

    def n_parameters(self):
        return 8


""" times = np.linspace(-100, 60, 1601)
p = [1, 1.2, 0.1195, 0.1003, 0.1089, 1.1357, 0.0453, 0]    # parameters
model = NfatModel()
values, markers = model.simulate(p, times)

print("markers = ", markers)

# Plot the results
plt.figure()
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.plot(times, values)
plt.show()

a = values[:, 0]
ac = values[:, 1]
anc = values[:, 2]
nc = values[:, 3]
an = p[0] - (a + ac + anc)
nn = p[1] - (nc + an + anc)
percentage_akap_bound_to_nfat = 100*(an+anc)/p[0]
percentage_nfat_on_membrane = 100*(an+anc)/p[1]

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
plt.show() """
