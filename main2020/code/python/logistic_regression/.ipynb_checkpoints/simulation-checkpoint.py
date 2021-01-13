import numpy as np
from scipy.special import comb
print(comb(5,2))
n = 50
omega = 0.5
xi = 0.5
k = 3

def make_likelihood_ratio_function(obs):
    def likelihood_ratio_function(omega, xi)
        out = 1
        for x in obs:
            out *= comb(k, x)* ( (1-omega) + 2**k * omega * (1-xi)**(k-x) * xi**x )
        return out
    return likelihood_ratio_function

def LRT_stat(likelihood_ratio_function):
    stat = 0
    for omega in range(101):
        for xi in range(1, 100):
            omega /= 100
            xi /= 100
            stat = max(likelihood_ratio_function(omega, xi))
    return stat

for _ in range(1000):