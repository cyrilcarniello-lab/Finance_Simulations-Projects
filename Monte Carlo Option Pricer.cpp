#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>


std::pair<double, double> Monte_Terminal(double s, double r, double phi, double t, std::mt19937& gen){
    std::normal_distribution<double> dist(0.0, 1.0);
    double z = dist(gen);
    double top = s*std::exp((r-(phi*phi)/2)*t+phi*std::sqrt(t)*z);
    double bot = s*std::exp((r-(phi*phi)/2)*t+phi*std::sqrt(t)*-z);
    return {top, bot};
}

std::pair<double, double> Monte_Pricer(double s, double k, double r, double phi, double t, double n, std::mt19937& gen){
    double sum = 0.0;
    double sqr_sum = 0.0;
    for (int i = 0; i<n; i++){
        auto result = Monte_Terminal(s, r, phi, t, gen);
        double top = result.first;
        double bot = result.second;
        double add = (std::max(bot-k, 0.0) + std::max(top-k, 0.0))/2;
        sum += add;
        sqr_sum += add*add;
    }
    double final_price = std::exp(-r*t)*(1/n)*sum;
    double stdt_var = (sqr_sum/n) - (sum/n)*(sum/n);
    double stdt_dev = std::sqrt(stdt_var);
    double stdt_err = stdt_dev/(std::sqrt(n))*std::exp(-r*t);
    
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