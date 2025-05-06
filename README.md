# SignalAnalysisTool

`SignalAnalysisTool` is a MATLAB class developed for analyzing the results of dynamic state and measurement estimation from Monte Carlo simulations. It is particularly tailored for power system applications such as synchronous generator models where the accuracy of state estimation under noise and uncertainty is critical.

## 🔍 Features

- Calculates **Mean Squared Error (MSE)** and **Mean Absolute Percentage Error (MAPE)** for each state and measurement
- Checks **convergence** of estimation results against a user-defined threshold
- Computes **statistical summaries** (mean ± 3σ) across Monte Carlo iterations
- Provides multiple **visualization tools** for estimation performance, including:
  - Time-series plots of estimated vs. true signals
  - Error bar plots across iterations
  - Error-vs-deviation plots
- Flexible option to **filter and analyze only convergent simulations**

## 📂 Structure

The class includes the following components:

### Properties

- `SimulationData` — stores original and estimated signals
- `ErrorStruct`, `MSEStruct`, `MAPE` — store error metrics
- `StatOfMonteCarlo` — contains statistical summaries of estimation results
- `Threshold`, `Deviations` — control convergence and noise/deviation analysis

### Public Methods

- `calculateMeanSquaredErrors()`
- `calculateMeanAbsolutePercentageError()`
- `calculateStatOfMonteCarlo()`
- `checkErrorConvergence(threshold)`
- `plotMonteCarloSimulationResults()`
- `plotErrorsBelowThreshold(threshold)`
- `plotMeanSquaredErrorTimesteps()`
- `plotMeanAbsolutePercentageErrorTimesteps()`
- And more...

## 🛠️ Usage

```matlab
% Load simulation results
simData = load('simulationData.mat');

% Create an instance
tool = SignalAnalysisTool(simData);

% Calculate error metrics
tool.calculateMeanSquaredErrors();
tool.calculateMeanAbsolutePercentageError();

% Check convergence with a custom threshold
[isX, isY] = tool.checkErrorConvergence(0.01);

% Visualize results
tool.plotMonteCarloSimulationResults('Threshold', true);
tool.plotErrorsBelowThreshold(0.05);
```

## 📊 Example Application

This tool is designed for post-processing results from state estimators such as:

* Extended Kalman Filter (EKF)
* Unscented Kalman Filter (UKF)
* Particle Filter (PF)

It can be integrated into power system studies, especially those involving:

* Generator model state estimation
* Robust estimation under noise and disturbances
* Benchmarking filtering performance across multiple trials

## 📎 Dependencies

* MATLAB R2020b or later
* Simulink (if data originates from Simulink models)
* Data structure format consistent with `x`, `x_hat`, `y`, `y_hat` fields

## 🚀 Getting Started

This repository is designed to work with MATLAB Simulink-based simulations for dynamic state estimation (e.g., synchronous generators). To use the `SignalAnalysisTool`, follow the steps below.

### 🔧 Requirements

- MATLAB R2020b or newer
- Simulink
- Estimation models for EKF, UKF, and PF (`proj_sim_ekf_wonoise`, `proj_sim_ukf_wonoise`, `proj_sim_pf_experiment`)
- Parameter script (e.g., `parameter_ekf_ukf_pf.m`)

---

## 🧪 Example Workflow (`main.m`)

A full Monte Carlo simulation and analysis can be done using the `main.m` script.

### Steps:
1. Select estimation method(s): EKF, UKF, PF.
2. Set up noise/deviation level for each iteration.
3. Run Simulink model for 100 Monte Carlo iterations.
4. Analyze estimation accuracy using:
   - `calculateMeanSquaredErrors()`
   - `calculateMeanAbsolutePercentageError()`
   - `checkErrorConvergence()`
   - `calculateStatOfMonteCarlo()`
5. Visualize:
   - Estimation vs. true states
   - Errors and convergence per iteration
   - Effect of deviation on MSE

### Run Simulation:

```matlab
% Inside main.m
simulationFileName = 'sim_pf_experiment';
parameter_ekf_ukf_pf; % Load system parameters

% Run simulation with selected filters (EKF, UKF, PF)
% Will execute Simulink model and analyze estimation accuracy
```
### Output:

* `.txt` file containing simulation log
* `.mat` file containing `SignalAnalysisTool` object per filter
* Multiple visual plots for diagnostic evaluation

## 📦 File Structure

```bash
.
├── main.m                    # Main simulation and analysis script
├── SignalAnalysisTool.m     # Core analysis class
├── parameter_ekf_ukf_pf.m   # System model and filter parameters
├── proj_sim_*.slx           # Simulink models for EKF, UKF, PF
├── *.mat                    # Output of simulations
└── *.txt                    # Logged simulation results
```

## 📘 License

This project is licensed under the MIT License — see the `LICENSE` file for details.

## 👨‍💻 Author

Azka Hariz Sartono<br>
Electrical Engineering — Universitas Gadjah Mada<br>
