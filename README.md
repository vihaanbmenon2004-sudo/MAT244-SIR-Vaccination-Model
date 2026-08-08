# MAT244 SIR Vaccination Model

This repository contains the MATLAB code used for our MAT244 final project on an extended Susceptible-Infected-Recovered (SIR) model with vaccination.

## Project Overview

The project studies how vaccination rate and vaccination timing affect the spread of an infectious disease. We compare:

- No vaccination
- Immediate constant vaccination
- Delayed vaccination
- Different combinations of transmission and vaccination rates
- Parameter sensitivity using ±20% changes in key model parameters

## Model

The extended SIR model includes a time-dependent vaccination rate that moves susceptible individuals directly into the immune population.

The baseline parameters used in the simulations are:

- Basic reproduction number: R0 = 2.5
- Recovery rate: gamma = 1/7 per day
- Transmission rate: beta = 2.5/7 per day
- Initial susceptible proportion: S(0) = 0.999
- Initial infected proportion: I(0) = 0.001
- Initial recovered proportion: R(0) = 0
- Simulation period: 180 days

## MATLAB Code

`sir_vaccination_project.m`

This script performs the numerical simulations using MATLAB's `ode45` solver.

It generates:

1. A comparison of immediate vaccination rates
2. A comparison of vaccination delays
3. A two-parameter sweep over transmission and vaccination rates
4. A ±20% parameter sensitivity analysis

## Figures

- `figure1_vaccination_rates.png` — Immediate vaccination-rate comparison
- `figure2_vaccination_delays.png` — Vaccination-delay comparison
- `figure3_heatmap.png` — Transmission-rate and vaccination-rate parameter sweep
- `figure4_sensitivity.png` — Parameter sensitivity analysis

## Running the Code

1. Open `sir_vaccination_project.m` in MATLAB.
2. Run the complete script.
3. Numerical results will appear in the MATLAB Command Window.
4. The four figures will be generated and saved automatically.
