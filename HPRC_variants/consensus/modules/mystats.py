from math import lgamma, gamma, exp, log, floor

from scipy.stats import poisson

# https://mathlandscape.com/poisson-distrib/

def continuous_poisson(x, lambda_):
    return exp(-lambda_) * lambda_**x / gamma(x+1)

def log_continuous_poisson(x, lambda_):
    if lambda_ == 0:
        # return log(0)
        return None
    # return -1 * lambda_ + x*log(lambda_) - log(gamma(x+1))
    return x*log(lambda_) - lambda_ - lgamma(x+1)


def log_continuous_cumulative_poisson(x, lambda_):
    # if x is not integer, return None
    if floor(x) != x:
        print("x is not integer")
        return None
    if lambda_ == 0:
        # return log(0)
        return None
    if x < 0:
        return 0
    if x > lambda_: #should use logsf
        return poisson.logsf(x, lambda_)
    elif x <= lambda_: #should use logcdf
        return poisson.logcdf(x, lambda_)



# test for x, lambda in 0 ~ 1

N = 1000
test_x = [ 1.0 / N * i for i in range(1, N+1)]
test_lambda = [ 1.0 / N * i for i in range(1, N+1)]
# get the minimum log continuous_poisson 
min_cp = min([log_continuous_poisson(x, lambda_) for x in test_x for lambda_ in test_lambda])

# heat map で表示
# x: lambda, y: x
# 0 ~ 1

# import matplotlib.pyplot as plt
# import numpy as np
# import seaborn as sns

# plt.figure()
# sns.heatmap(np.array([[log_continuous_poisson(x, lambda_)  for x in test_x] for lambda_ in test_lambda]))
# plt.xlabel("lambda")
# plt.ylabel("x")
# plt.title("log continuous poisson")
# # colorbar
# plt.savefig("./heatmap.png")
