import numpy as np
from scipy.special import comb
from scipy.stats import binom, chi2, beta, norm, multivariate_normal
from scipy import optimize
import matplotlib.pyplot as plt
import math

class binomial_simulator():
    def __init__(self, n, k, omega, xi, numerator_power = 2/3, denominator_power = 1/3, sampling_N = 100, M = 10000, M_gaussian = 20000):
        self.n = n
        self.k = k
        self.omega = omega
        self.xi = xi

        self.numerator_power = numerator_power
        self.denominator_power = denominator_power

        # sampling number for computing the LRT and the integrated likelihood
        self.sampling_N = sampling_N

        # repeat times
        self.M = M
        # repeat times for guassian process simulation
        self.M_gaussian = M_gaussian

        print("n", self.n)
        print("k", self.k)
        print("omega", self.omega)
        print("xi", self.xi)

    def generate_data(self):
        # data generation
        tmp1 = binom.rvs(self.k, 0.5, size = self.n)
        tmp2 = binom.rvs(self.k, self.xi, size = self.n)
        switch = binom.rvs(1, self.omega, size = self.n)
        self.obs = switch * tmp2 + (1-switch) * tmp1

        self.dense_obs = [0 for _ in range(self.k+1)]
        for i in self.obs:
            self.dense_obs[i] += 1

    def log_likelihood_ratio_function(self, omega, xi):

        if omega < 0 or omega >1 or xi <= 0 or xi >= 1:
            return -1e6
        out = 0
        for i in range(self.k + 1):
            tmp = self.dense_obs[i] * math.log( (1-omega) + 2**self.k * omega * (1-xi)**(self.k-i) * xi**i )
            out += tmp
        return out

    def do_one_simul(self):
        log_LRT = -1e6
        gFBF_numerator = 0
        gFBF_denominator = 0
        omega_list = [i/self.sampling_N for i in range(self.sampling_N + 1 )]
        xi_list = [i/self.sampling_N for i in range(1, self.sampling_N)]
        for omega in omega_list:
            for xi in xi_list:
                tmp = self.log_likelihood_ratio_function(omega, xi)
                # LRT
                if tmp > log_LRT:
                    log_LRT = tmp

                # gFBF
                # prior for oemga is beta(2,2), prior for xi is uniform(0,1)
                prior = 6 * omega * (1-omega)
                gFBF_numerator += math.exp(tmp * self.numerator_power) * prior
                gFBF_denominator += math.exp(tmp * self.denominator_power) * prior

        log_gFBF = math.log( gFBF_numerator / gFBF_denominator )

        out = {}
        out["LRT"] = 2 * log_LRT
        out["gFBF"] = 2 / (self.numerator_power - self.denominator_power) * log_gFBF + 1/ ( self.numerator_power - self.denominator_power) * math.log(self.numerator_power / self.denominator_power)

        return out["LRT"], out["gFBF"]

    def gaussian_process_quantile(self):
        the_f = [0 for _ in range(self.k+1)]
        for i in range(self.k+1):
            the_f[i] = 2**(-self.k) * comb(self.k, i)
        the_Sigma = np.zeros((self.k+1, self.k+1))
        for i in range(self.k+1):
            for j in range(self.k+1):
                the_Sigma[i,j] = - the_f[i] * the_f[j]
                if i == j:
                    the_Sigma[i,j] += the_f[i]

        def P_kx(k, x, phi):
            if phi == 0:
                return 2 * x - k
            return ( (1+phi)**x * (1-phi)**(k-x) - 1) / phi

        # simulate the gaussian process to approximate the distribution of -2 log LRT
        def simu_gauss_process():
            data = multivariate_normal.rvs(mean = None, cov = the_Sigma, size = 1)
            phi_list = [ i / 1000 * 2 - 1 for i in range(1001)]
            T_squared = 0
            for phi in phi_list:
                S = 0
                for j in range(self.k+1):
                    S += data[j] * P_kx(self.k, j, phi)
                if phi < 0 :
                    S *= -1
                S = max(S, 0)
                tmp = S ** 2 / P_kx(self.k, self.k, phi**2)
                if tmp > T_squared:
                    T_squared = tmp
            return T_squared
        gauss_process_stats = [0 for _ in range( self.M_gaussian)]
        for i in range(self.M_gaussian):
            gauss_process_stats[i] = simu_gauss_process()
        self.gaussian_process_critical_value = np.quantile(gauss_process_stats, 0.95)

        #print("gauss process quantile:")
        #for q in [0.25, 0.5, 0.75, 0.90, 0.95, 0.99]:
        #    print(q, np.quantile(gauss_process_stats, q))


    def do_it(self):

        LRT_stats = [0 for _ in range (self.M)]
        gFBF_stats = [0 for _ in range (self.M)]

        LRT_rejects = [0 for _ in range (self.M)]
        gFBF_rejects = [0 for _ in range (self.M)]

        gFBF_critical_value = chi2.ppf(0.95, df = 1)
        self.gaussian_process_quantile()

        for i in range(self.M):
            self.generate_data()
            LRT_stats[i], gFBF_stats[i] = self.do_one_simul()
            if LRT_stats[i] > self.gaussian_process_critical_value:
                LRT_rejects[i] = 1
            if gFBF_stats[i] > gFBF_critical_value:
                gFBF_rejects[i] = 1

        #print("chi-squared quantile:")
        #for q in [0.25, 0.5, 0.75, 0.90, 0.95, 0.99]:
        #    print(q, chi2.ppf(q, df = 1) )

        #print("LRT quantile:")
        #for q in [0.25, 0.5, 0.75, 0.90, 0.95, 0.99]:
        #    print(q, np.quantile(LRT_stats, q))
        #print("fFBF quantile:")
        #for q in [0.25, 0.5, 0.75, 0.90, 0.95, 0.99]:
        #    print(q, np.quantile(gFBF_stats, q))

        print("numerator_power", self.numerator_power)
        print("denominator_power", self.denominator_power)

        print("LRT_power", np.mean(LRT_rejects))
        print("gFBF_power", np.mean(gFBF_rejects))


if __name__ == "__main__" :
    import time
    now_time = time.time()
    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.5, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()
    print("time", time.time()-now_time)

    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.6, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.7, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.8, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.9, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.25, xi = 0.6, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.25, xi = 0.7, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.25, xi = 0.8, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.25, xi = 0.9, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.75, xi = 0.6, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.75, xi = 0.7, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.75, xi = 0.8, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.75, xi = 0.9, numerator_power = 4/5, denominator_power = 1/5)
    simulator.do_it()






    now_time = time.time()
    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.5, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()
    print("time", time.time()-now_time)

    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.6, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.7, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.8, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.5, xi = 0.9, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.25, xi = 0.6, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.25, xi = 0.7, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.25, xi = 0.8, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.25, xi = 0.9, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.75, xi = 0.6, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.75, xi = 0.7, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.75, xi = 0.8, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

    simulator = binomial_simulator (n=50, k = 4, omega = 0.75, xi = 0.9, numerator_power = 3/5, denominator_power = 2/5)
    simulator.do_it()

