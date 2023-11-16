from math import gamma, exp, log

def continuous_poisson(x, lambda_):
    return exp(-lambda_) * lambda_**x / gamma(x+1)

def log_continuous_poisson(x, lambda_):
    if lambda_ == 0:
        # return log(0)
        return None
    return -1 * lambda_ + x*log(lambda_) - log(gamma(x+1))