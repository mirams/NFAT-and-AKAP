import numpy as np
import pints
from scipy import interpolate
from model import NfatModel

# For debug only
# import matplotlib.pyplot as plt

akap_nfatp_times = np.genfromtxt(
    'data/NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(0))
akap_nfatp_means = np.genfromtxt(
    'data/NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(2))

# Scale down such that the mean of the first 4 points is at 100%.
akap_nfatp_means = akap_nfatp_means/np.mean(akap_nfatp_means[0:3])
# akap_nfatp_sems = np.genfromtxt(
#    'data/NFAT-AKAP disasociation-F7E.csv', delimiter=',', skip_header=1, usecols=(3))


class SpecialLogLikelihood(pints.LogPDF):

    def __call__(self, params):
        assert(len(params) == 5)

        for param in params:
            if param < 0:
                print("Less than 0")
            if param > 2:
                print("More than 2!")

        times = np.linspace(-500, 60, 56001)
        model = NfatModel()
        p = np.array([1, 2])

        values, markers = model.simulate(
            np.concatenate([p, params, [0.1, 0.0]]), times)

        if len(values) != len(times):
            log_likelihood = -1e9*sum(params**2)
            return

        # correct_answers = [24, 11, 26, 17, 8, 12] # These were the original confocal estimates
        # These are the new 3D scan estimates for NFAT locations
        correct_answers = [12, 5.6, 26, 17, 8, 12]
        sigmas = [1.7, 1.7, 9, 6, 3, 3]

        log_likelihood = 0
        for i in range(0, len(correct_answers)):
            log_likelihood = log_likelihood - \
                (markers[i] - correct_answers[i]) ** 2
        log_likelihood = log_likelihood * (1.0/(2.0*(sigmas[i]*sigmas[i])))
        # print('LL = ', log_likelihood)

        a = values[:, 0]
        ac = values[:, 1]
        anc = values[:, 2]
        nc = values[:, 3]
        an = p[0] - (a + ac + anc)
        nn = p[1] - (nc + an + anc)
        total_akap_with_nfat = an+anc
        proportion = total_akap_with_nfat/np.max(total_akap_with_nfat)

        # For interpolating simulation results to the data times
        akap_nfat_prediction_interpolator = interpolate.interp1d(
            times, proportion)

        # Add a constraint for nuclear proportion of NFAT at the end to be about 80%
        # Weight on fitting this (smaller = more weight)
        nuc_prop_sigma = 1.0
        nuclear_proportion = nn/p[1]
        min_20_time_index = np.where(np.abs(times-20) < 1e-12)
        nuclear_proportion_time_20 = nuclear_proportion[min_20_time_index[0]][0]
        log_likelihood = log_likelihood - 0.5*np.log(2*np.pi*nuc_prop_sigma**2)-(
            1.0/(2.0*nuc_prop_sigma**2))*(nuclear_proportion_time_20 - 0.80)**2

        # std_of_timecourse_data = len(
        #    akap_nfatp_times)*(akap_nfatp_sems[i]*akap_nfatp_sems[i])
        # print('Sigma timecourse = ', std_of_timecourse_data)

        timecourse_sigma = 0.07
        sum_sq = 0
        for i in range(1, len(akap_nfatp_times)):
            sum_sq = sum_sq+(akap_nfat_prediction_interpolator(
                akap_nfatp_times[i]) - akap_nfatp_means[i]) ** 2
        log_likelihood = log_likelihood - \
            (len(akap_nfatp_times)/2)*np.log(2*np.pi *
                                             timecourse_sigma**2) - 1.0 / (2.0*timecourse_sigma**2) * sum_sq
        # print('t=', akap_nfatp_times[i], ' Prediction = ', akap_nfat_prediction_interpolator(
        #    akap_nfatp_times[i]), ' Data = ', akap_nfatp_means[i], ' LL = ', log_likelihood)

        # plt.figure('From within Special Log Likelihood')
        # plt.xlabel('Time')
        # plt.ylabel('Proportion of original AKAP:NFATp co-localised (%)')
        # plt.plot(akap_nfatp_times, 100*akap_nfatp_means, 'x')
        # plt.plot(times, 100*(an+anc)/np.max(an+anc))
        # plt.show()

        # Add a penalty for k2 being too small and AC never dissociating back to A (and ANC never dissociating back to AN)
        penalty_sigma = 0.1
        neg_time_index = np.where(np.abs(times) < 1e-12)
        target_end_an = an[neg_time_index[0]][0]
        log_likelihood = log_likelihood - 1.0 / \
            (2.0*penalty_sigma**2)*(target_end_an - an[-1])**2

        return log_likelihood

    def n_parameters(self):
        return 5
