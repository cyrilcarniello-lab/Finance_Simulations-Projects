import math
def black_pricer(s, k, r, phi, t):
    d2=(math.log(s/k)+(r-(phi**2)/2)*t)/(phi*math.sqrt(t))
    d1=d2+phi*math.sqrt(t)
    return s*N(d1) - k*math.exp(-r*t)*N(d2)
    


def N(x):
    return 0.5* (1 + math.erf(x/math.sqrt(2)))