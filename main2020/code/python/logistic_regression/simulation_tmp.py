import numpy as np
from scipy.special import comb
from scipy.stats import binom, chi2, beta, norm, multivariate_normal, cauchy
from scipy import optimize
import matplotlib.pyplot as plt
import math
import statsmodels.api as sm
import random

class logistic_simulator():
    def __init__(self, n, p, p0, true_nu, true_xi, numerator_power = 2/3, denominator_power = 1/3, sampling_N = 1000, M = 10000):

        # the null hypothesis is xi = 0

        self.n = n
        self.p = p
        self.p0 = p0
        self.true_nu = np.array(true_nu)
        self.true_xi = np.array(true_xi)
        self.true_beta = np.concatenate((self.true_nu, self.true_xi), axis = 0)

        self.numerator_power = numerator_power
        self.denominator_power = denominator_power

        # sampling number for computing the integrated likelihood
        self.sampling_N = sampling_N

        # repeat times
        self.M = M

        print("n", self.n)
        print("true_beta", self.true_beta)
    def sigmoid(self, x):
        return 1/(1+math.exp(-x))

    def generate_data(self):
        # data generation
        X1 = norm.rvs(0, 1, size = self.n * self.p0)
        X1 = X1.reshape((self.n, self.p0))
        X2 = norm.rvs(0, 1, size = self.n * (self.p - self.p0))
        X2 = X2.reshape((self.n, self.p - self.p0))
        self.X1 = X1
        self.X = np.concatenate((X1, X2), axis = 1)

        tmp = self.X.dot(self.true_beta)
        self.y = [ int(random.random() < self.sigmoid(i)) for i in tmp]

    def do_one_simul(self):
        out = {}
        gFBF_full_numerator = 0
        gFBF_full_denominator = 0
        gFBF_null_numerator = 0
        gFBF_null_denominator = 0

        sm_model_full = sm.Logit(self.y, self.X)
        try:
            #sm_result_full = sm_model_full.fit(method="newton", callback = self.my_callback)
            #sm_result_full = sm_model_full.fit(method="newton")
            full_llf = sm_model_full.fit(method="bfgs").llf
        except:
            full_error = True
            full_llf = 0
        else:
            full_error = False

        sm_model_null = sm.Logit(self.y, self.X1)
        try:
            #sm_result_null = sm_model_null.fit(method="newton", callback = self.my_callback)
            #sm_result_null = sm_model_null.fit(method="newton")
            null_llf = sm_model_null.fit(method="bfgs").llf
        except:
            null_error = True
            null_llf = 0
        else:
            null_error = False

        out["LRT"] = 2 * (full_llf - null_llf)

        params_total_full = cauchy.rvs(0, 2.5, size = self.p * self.sampling_N)
        params_total_null = cauchy.rvs(0, 2.5, size = self.p0 * self.sampling_N)

        for i in range(self.sampling_N):
            # gFBF
            #params = cauchy.rvs(0, 2.5, size = self.p)

            params_full = params_total_full[(i*self.p): ((i+1)*self.p)]
            loglikelihood_ratio = sm_model_full.loglike(params_full) + self.n * math.log(2)
            gFBF_full_numerator += math.exp(loglikelihood_ratio * self.numerator_power)
            gFBF_full_denominator += math.exp(loglikelihood_ratio * self.denominator_power)

            params_null = params_total_null[(i*self.p0): ((i+1)*self.p0)]
            loglikelihood_ratio = sm_model_null.loglike(params_null) + self.n * math.log(2)
            gFBF_null_numerator += math.exp(loglikelihood_ratio * self.numerator_power)
            gFBF_null_denominator += math.exp(loglikelihood_ratio * self.denominator_power)

        log_gFBF = math.log( (gFBF_full_numerator / gFBF_full_denominator )/(gFBF_null_numerator / gFBF_null_denominator ) )

        out["gFBF"] = 2 / (self.numerator_power - self.denominator_power) * log_gFBF + (self.p - self.p0)/ ( self.numerator_power - self.denominator_power) * math.log(self.numerator_power / self.denominator_power)

        return out["LRT"], out["gFBF"]

    def do_it(self):

        LRT_stats = [0 for _ in range (self.M)]
        gFBF_stats = [0 for _ in range (self.M)]

        LRT_rejects = [0 for _ in range (self.M)]
        gFBF_rejects = [0 for _ in range (self.M)]

        gFBF_critical_value = chi2.ppf(0.95, df = self.p - self.p0)
        LRT_critical_value = chi2.ppf(0.95, df = self.p - self.p0)

        for i in range(self.M):
            self.generate_data()
            LRT_stats[i], gFBF_stats[i] = self.do_one_simul()
            if LRT_stats[i] > LRT_critical_value:
                LRT_rejects[i] = 1
            if gFBF_stats[i] > gFBF_critical_value:
                gFBF_rejects[i] = 1

        print("chi-squared quantile:")
        for q in [0.25, 0.5, 0.75, 0.90, 0.95, 0.99]:
            print(q, chi2.ppf(q, df = 1) )

        print("LRT quantile:")
        for q in [0.25, 0.5, 0.75, 0.90, 0.95, 0.99]:
            print(q, np.quantile(LRT_stats, q))
        print("fFBF quantile:")
        for q in [0.25, 0.5, 0.75, 0.90, 0.95, 0.99]:
            print(q, np.quantile(gFBF_stats, q))

        print("LRT_power", np.mean(LRT_rejects))
        print("gFBF_power", np.mean(gFBF_rejects))
        return {"LRT": np.mean(LRT_rejects), "gFBF": np.mean(gFBF_rejects)}

if __name__ == "__main__" :
    import time
    now_time = time.time()
    simulator = logistic_simulator(n =50, p=3, p0=2, true_nu=[1, 1], true_xi = [10] * 1)
    simulator.do_it()
    print("time", time.time()-now_time)


    #simulator = normal_simulator (n=50, omega = 0.5, xi = 0.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.5, xi = 1)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.5, xi = 1.5)
    #simulator.do_it()

    #simulator = normal_simulator (n=50, omega = 0.5, xi = 2)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.5, xi = 0)
    #simulator.do_it()

    #simulator = normal_simulator (n=100, omega = 0.5, xi = 0.5)
    #simulator.do_it()
