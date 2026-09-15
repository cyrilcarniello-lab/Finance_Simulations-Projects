[Monte_Carlo_EU_Option.md](https://github.com/user-attachments/files/32250581/Monte_Carlo_EU_Option.md)
# Monte Carlo European Option Pricer

A European call option pricer built from first principles in both **Python** and **C++**, using Monte Carlo simulation under the Black-Scholes / Geometric Brownian Motion model. Every formula used — including the risk-neutral pricing framework, the Itô-corrected GBM terminal price formula, and the Black-Scholes closed form — was derived by hand rather than taken as given, and validated against the exact closed-form solution.

## Overview

The project prices a vanilla European call option, `max(S_T - K, 0)`, by simulating many possible future stock prices under risk-neutral Geometric Brownian Motion and averaging the resulting discounted payoffs:

```
Price ≈ e^(-rT) · (1/N) · Σ max(S_T,i - K, 0)
```

where each simulated terminal price is generated from:

```
S_T = S0 · exp((r - σ²/2)·T + σ·√T·Z),    Z ~ N(0,1)
```

The `-σ²/2` term is the Itô correction, arising from applying Itô's lemma to the geometric Brownian motion SDE `dS = rS·dt + σS·dW` — derived explicitly rather than assumed.

## Features

- **Vanilla Monte Carlo pricer** — simulates terminal stock prices and averages discounted payoffs.
- **Custom random normal generation** — standard normal draws implemented via the Box-Muller transform (Python) and `std::mt19937` + `std::normal_distribution` (C++), properly seeded with `std::random_device` to avoid the classic fixed/low-resolution seeding mistake.
- **Standard error estimation** — tracks the Monte Carlo estimator's own uncertainty (`σ_payoff/√N`, discounted), not just the point estimate.
- **Closed-form Black-Scholes validator** — an independently implemented exact pricer (`N(x)` computed via the error function, `erf`) used to check the simulation's accuracy.
- **Convergence analysis** — price vs. `N` with confidence bands, and a log-log plot of absolute error vs. `N`, empirically confirming the theoretical `O(1/√N)` convergence rate (slope ≈ -1/2).
- **Variance reduction (antithetic variates)** — pairs each draw `Z` with its mirror `-Z` and averages the resulting payoffs, reducing the standard error beyond what additional independent sampling alone would achieve.
- **Ported to C++** — same model and variance-reduction technique reimplemented in C++, for a direct language/performance comparison.

## Results

For `S0 = 100, K = 100, r = 0.05, σ = 0.2, T = 1`:

| Method | Price | Notes |
|---|---|---|
| Black-Scholes (exact) | ≈ $10.45 | Closed-form reference |
| Monte Carlo (vanilla) | ≈ $10.4x | N = 100,000, standard error ≈ $0.046 |
| Monte Carlo (antithetic variates) | ≈ $10.44 | N = 100,000, standard error ≈ $0.023 — roughly half, from genuine negative correlation between paired draws, not just added samples |

The convergence plot confirms the Monte Carlo error shrinks at the theoretically predicted `1/√N` rate: the log-log slope of absolute error vs. `N` matches `-1/2`.

## Project structure

```
.
├── black_scholes_pricer.py       # Closed-form Black-Scholes validator (Python)
├── monte_carlo_pricer.py         # Monte Carlo pricer, convergence plots, variance reduction (Python)
├── monte_carlo_pricer.cpp        # C++ port of the Monte Carlo pricer
└── README.md
```

## Running it

**Python:**
```bash
python3 monte_carlo_pricer.py
```

**C++:**
```bash
clang++ monte_carlo_pricer.cpp -o pricer
./pricer
```

## Background

This project was built as a learning exercise covering: risk-neutral option pricing, Geometric Brownian Motion and Itô's lemma, the Box-Muller transform for normal random variate generation, Monte Carlo convergence theory, variance reduction techniques, and C++ fundamentals (`<random>`, references, `std::pair`) coming from a C background.
