import numpy as np
from scipy.special import comb
from scipy.stats import binom, chi2, beta, norm, multivariate_normal
from scipy import optimize
import matplotlib.pyplot as plt
import math

class normal_simulator():
    def __init__(self, n, omega, xi, numerator_power = 2/3, denominator_power = 1/3, sampling_N = 100, M = 10000, M_gaussian = 20000):
        self.n = n
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

        self.LRT_critical_value = math.log(math.log(n)) - math.log(2* math.pi) - 2 * math.log(- math.log(0.95))

        print("n", self.n)
        print("omega", self.omega)
        print("xi", self.xi)

    def generate_data(self):
        # data generation
        tmp1 = norm.rvs(0, 1, size = self.n)
        tmp2 = norm.rvs(self.xi, 1, size = self.n)
        switch = binom.rvs(1, self.omega, size = self.n)
        self.obs = switch * tmp2 + (1-switch) * tmp1

    def log_likelihood_ratio_function(self, omega, xi):

        if omega < 0 or omega >1:
            return -1e6
        out = 0
        for i in range(self.n):
            tmp =  math.log( (1-omega) + omega * math.exp( self.obs[i] * xi - xi * xi /2 ) )
            out += tmp
        return out

    def do_one_simul(self):
        log_LRT = -1e6
        gFBF_numerator = 0
        gFBF_denominator = 0
        omega_list = [i/self.sampling_N for i in range(self.sampling_N + 1 )]
        xi_list = [(i/self.sampling_N - 1/2) * 5 for i in range(1, self.sampling_N)]
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


    def do_it(self):

        LRT_stats = [0 for _ in range (self.M)]
        gFBF_stats = [0 for _ in range (self.M)]

        LRT_rejects = [0 for _ in range (self.M)]
        gFBF_rejects = [0 for _ in range (self.M)]

        gFBF_critical_value = chi2.ppf(0.95, df = 1)

        for i in range(self.M):
            self.generate_data()
            LRT_stats[i], gFBF_stats[i] = self.do_one_simul()
            if LRT_stats[i] > self.LRT_critical_value:
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

        print("LRT_power", np.mean(LRT_rejects))
        print("gFBF_power", np.mean(gFBF_rejects))


if __name__ == "__main__" :
    import time
    now_time = time.time()
    #simulator = normal_simulator (n=50, omega = 0.5, xi = 0)
    #simulator.do_it()
    print("time", time.time()-now_time)


    #simulator = normal_simulator (n=50, omega = 0.5, xi = 0.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.5, xi = 1)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.5, xi = 1.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.5, xi = 2)
    #simulator.do_it()




    #simulator = normal_simulator (n=50, omega = 0.25, xi = 0.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.25, xi = 1)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.25, xi = 1.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.25, xi = 2)
    #simulator.do_it()




    #simulator = normal_simulator (n=50, omega = 0.75, xi = 0.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.75, xi = 1)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.75, xi = 1.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.75, xi = 2)
    #simulator.do_it()





    simulator = normal_simulator (n=100, omega = 0.5, xi = 0)
    simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.5, xi = 0.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.5, xi = 1)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.5, xi = 1.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.5, xi = 2)
    #simulator.do_it()




    #simulator = normal_simulator (n=100, omega = 0.25, xi = 0.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.25, xi = 1)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.25, xi = 1.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.25, xi = 2)
    #simulator.do_it()




    #simulator = normal_simulator (n=100, omega = 0.75, xi = 0.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.75, xi = 1)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.75, xi = 1.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.75, xi = 2)
    #simulator.do_it()
