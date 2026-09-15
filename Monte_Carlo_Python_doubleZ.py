import math
from random import random
import matplotlib.pyplot as plt
from Black_Scholes_Pricer import black_pricer

def Monte_pricer(s, k, r, phi, t, n):
    
    final_sum = 0
    val = 0
    m = n
    stdt_coun = 0
    while m != 0:
        pos, neg = Terminal_val(s, r, phi, t)
        val = (max(pos-k, 0) + max(neg-k, 0))/2
        final_sum += val
        stdt_coun += val**2
        m -= 1

    final_price = math.exp(-r*t)*(1/n)*final_sum
    stdt_var = (stdt_coun/n) - (final_sum/n)**2
    stdt_dev = math.sqrt(stdt_var)
    stdt_err = stdt_dev/(math.sqrt(n))*math.exp(-r*t)
    #print(stdt_var)
    #print(stdt_dev)
    return final_price, stdt_err


def Terminal_val(s, r, phi, t):
    A = random()
    B = random()
    Z = math.sqrt(-2*math.log(A))*math.cos(2*math.pi*B)
    return s*math.exp((r-(phi**2)/2)*t+phi*math.sqrt(t)*Z), s*math.exp((r-(phi**2)/2)*t+phi*math.sqrt(t)*-Z)

print(Monte_pricer(100, 100, 0.05, 0.2, 1, 100000))



N_values = []
Err_values_top = []
Err_values_bot = []
Err_time = []
prices = []
j = 1

exact = black_pricer(100, 100, 0.05, 0.2, 1)

while j < 200000:
    j *= 5
    N_values.append(j)                      #creates list with no of sim vals up to 200000; logarithmic, to reduce runtime

for i in N_values:                          #lists to later insert into the plot(s)
    price, stdt_error = Monte_pricer(100, 100, 0.05, 0.2, 1, i)
    prices.append(price)
    Err_values_top.append(price + 2*stdt_error)
    Err_values_bot.append(price - 2*stdt_error)
    Err_time.append(abs(price-exact))       

plt.plot(N_values, prices)                  #creating plot to compare the Terminal_vals based on num of simulations
plt.plot(N_values, Err_values_top)          #shows the upper boundary based on the calculated stadart error   
plt.plot(N_values, Err_values_bot)          #shows same as above, but for the lower boundary
plt.xscale('log')
plt.xlabel("No. of Simulations")
plt.ylabel("Cummulated Price")
plt.show()                                  #shows the created plot


plt.plot(N_values, Err_time)                #creating plot to show how the difference betweeen the exact val(through black-Scholes) adn the Monte_Carlo val gets smaller with no of simulations
plt.xscale('log')
plt.yscale('log')
plt.xlabel("No. of Simulations")
plt.ylabel("Error")
plt.show()