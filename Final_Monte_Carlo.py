import math
from random import random
import matplotlib.pyplot as plt

def Monte_pricer(s, k, r, phi, t, n):
    
    final_sum = 0                               #final sum of all simulated values
    val = 0
    m = n
    stdt_coun = 0                               #final sum of all squared values(needed for stdt_error calculation)
    while m != 0:                               #loop running through the simulations
        pos, neg = Terminal_val(s, r, phi, t)   #getting the return values through the Terminal_val helper function -> both vals with z and -z calculated
        val = (max(pos-k, 0) + max(neg-k, 0))/2 #getting th avg of the two return vals
        final_sum += val
        stdt_coun += val**2                     
        m -= 1

    final_price = math.exp(-r*t)*(1/n)*final_sum    #calculation of final price through the GBM formula
    stdt_var = (stdt_coun/n) - (final_sum/n)**2
    stdt_dev = math.sqrt(stdt_var)
    stdt_err = stdt_dev/(math.sqrt(n))*math.exp(-r*t)
    return final_price, stdt_err                    #returns both standart error && final price 


def Terminal_val(s, r, phi, t):                     #helper function, to determine both return values of each simulation run
    A = random()                                    #random vals betweeen 0 && 1  
    B = random()                                    
    Z = math.sqrt(-2*math.log(A))*math.cos(2*math.pi*B)     #simulation of a z in a N(0,1)->based on Box-Mueller
    return s*math.exp((r-(phi**2)/2)*t+phi*math.sqrt(t)*Z), s*math.exp((r-(phi**2)/2)*t+phi*math.sqrt(t)*-Z)   #returning both values of the positive and negative z -> halves stdt_err 

print(Monte_pricer(100, 100, 0.05, 0.2, 1, 100000))

def black_pricer(s, k, r, phi, t):          #black scholes pricer; to calculate the vals for the Err_time Model
    d2=(math.log(s/k)+(r-(phi**2)/2)*t)/(phi*math.sqrt(t))
    d1=d2+phi*math.sqrt(t)
    return s*N(d1) - k*math.exp(-r*t)*N(d2)
    


def N(x):                                   #helper function to the black scholes func above
    return 0.5* (1 + math.erf(x/math.sqrt(2)))


N_values = []                                       #plots based on values from above
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



