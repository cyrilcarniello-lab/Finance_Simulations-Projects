#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>


std::pair<double, double> Monte_Terminal(double s, double r, double phi, double t, std::mt19937& gen){              //Monte_Terminal evaluates the terminal values of the option with a randomized z/-z value
    std::normal_distribution<double> dist(0.0, 1.0);                                //marsenne twister allows us to use the already implemented normal distribution function -> no direct box-mueller based calculation needed
    double z = dist(gen);
    double top = s*std::exp((r-(phi*phi)/2)*t+phi*std::sqrt(t)*z);
    double bot = s*std::exp((r-(phi*phi)/2)*t+phi*std::sqrt(t)*-z);
    return {top, bot};
}

std::pair<double, double> Monte_Pricer(double s, double k, double r, double phi, double t, double n, std::mt19937& gen){
    double sum = 0.0;                                                               //final sum of all simulated values
    double sqr_sum = 0.0;
    for (int i = 0; i<n; i++){                                                      //loop running through the Monte_Terminal-simulations n-times
        auto result = Monte_Terminal(s, r, phi, t, gen);                            //Monte_Terminal returns both prices for the simulated z and -z values
        double top = result.first;
        double bot = result.second;
        double add = (std::max(bot-k, 0.0) + std::max(top-k, 0.0))/2;               //double add takes the average of the values top bot; after determining if they even turn profit; else their counted as 0
        sum += add;
        sqr_sum += add*add;                                                         //the sqr_sum sums all squared add values; to later determine the standart error
    }
    double final_price = std::exp(-r*t)*(1/n)*sum;                                  //final price calculation gets the average from the sum and at the same time continously compounding the rate to the final value
    double stdt_var = (sqr_sum/n) - (sum/n)*(sum/n);
    double stdt_dev = std::sqrt(stdt_var);
    double stdt_err = stdt_dev/(std::sqrt(n))*std::exp(-r*t);                       //stdt error calculation based on earlier calculated standart variance and deviation
    
    return {final_price, stdt_err};
}


int main() {
    std::random_device rd;
    std::mt19937 gen(rd());                             //mersenne twister
    std::normal_distribution<double> dist(0.0, 1.0);
    double z = dist(gen);
    double s = 100;
    double k = 100;
    double r = 0.05;
    double phi = 0.2;
    double t = 1;
    double n = 100000;
    auto result = Monte_Pricer(s, k, r, phi, t, n, gen);
    std::cout << result.first << std::endl;
    std::cout << result.second << std::endl;
    return 0;
}