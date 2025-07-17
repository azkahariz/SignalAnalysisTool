# SignalAnalysisTool

`SignalAnalysisTool` is a MATLAB class developed for the post-processing and analysis of dynamic state and measurement estimation results obtained from Monte Carlo simulations. This tool is specifically designed for power system applications, such as synchronous generator models, where accuracy under varying noise conditions and uncertainties is critical.

## 🔍 Features

* Computes **Mean Squared Error (MSE)** and **Mean Absolute Percentage Error (MAPE)** for each state variable and measurement.
* Evaluates **convergence** of estimation results based on a user-defined error threshold.
* Provides **statistical summaries** (mean ± 3σ) across Monte Carlo iterations.
* Offers comprehensive **visualization tools**, including:

  * Time-series plots comparing estimated and true signals.
  * Error bar plots across different iterations.
  * Error versus deviation plots.
  * Comparative performance plots across various noise scenarios.
* Enables **selective analysis** of only convergent simulation results.

## 📂 Class Structure

### Properties

* `SimulationData` — Stores true and estimated state and measurement signals.
* `ErrorStruct` and `MSEStruct` — Store error data and MSE calculations.
* `MAPE` — Stores Mean Absolute Percentage Error results.
* `StatOfMonteCarlo` — Contains statistical summaries (mean and standard deviation).
* `Threshold` — Convergence threshold definition.
* `Deviations` — Stores deviation data for error analysis.

### Key Public Methods

* `calculateMeanSquaredErrors()` — Computes MSE.
* `calculateMeanAbsolutePercentageError()` — Computes MAPE.
* `calculateStatOfMonteCarlo()` — Calculates mean and standard deviation of results.
* `checkErrorConvergence(threshold)` — Determines convergence status.
* `plotMonteCarloSimulationResults()` — Plots time-series estimation results.
* `plotErrorsBelowThreshold(threshold)` — Visualizes errors within threshold.
* `plotMeanSquaredErrorTimesteps()` — Visualizes MSE over simulation time.
* `plotMeanAbsolutePercentageErrorTimesteps()` — Visualizes MAPE over simulation time.

## 🛠️ Basic Usage Example

```matlab
% Load simulation results
simData = load('simulationData.mat');

% Instantiate the analysis tool
tool = SignalAnalysisTool(simData);

% Perform error analysis
tool.calculateMeanSquaredErrors();
tool.calculateMeanAbsolutePercentageError();

% Evaluate convergence
[isX, isY] = tool.checkErrorConvergence(0.01);

% Visualize estimation results
tool.plotMonteCarloSimulationResults('Threshold', true);
tool.plotErrorsBelowThreshold(0.05);
```

## 📊 Application Scope

This tool is intended for analyzing post-simulation results from state estimation methods such as:

* Extended Kalman Filter (EKF)
* Unscented Kalman Filter (UKF)
* Particle Filter (PF)

It supports power system studies focused on:

* Synchronous generator state estimation
* Performance evaluation under various noise disturbances
* Benchmarking filtering methods in dynamic systems

## 📎 Dependencies

* MATLAB R2020b or later
* Simulink (if using Simulink-based simulation data)
* Simulation data formatted with required fields (`x`, `x_hat`, `y`, `y_hat`)

## 🚀 Getting Started

This repository is designed for use with Simulink-based dynamic system simulations (e.g., synchronous generators). Follow the steps below to utilize the `SignalAnalysisTool` class:

### Requirements

* MATLAB R2020b or newer
* Simulink
* Predefined estimation models: `proj_sim_ekf_wonoise`, `proj_sim_ukf_wonoise`, `proj_sim_pf_experiment`
* Parameter initialization script (e.g., `parameter_ekf_ukf_pf.m`)

## 🧪 Example Workflow (from `main.m` script)

### Process Overview

1. Select estimation methods: EKF, UKF, or PF.
2. Define noise or deviation parameters.
3. Execute Simulink model for Monte Carlo iterations.
4. Analyze estimation results via:

   * `calculateMeanSquaredErrors()`
   * `calculateMeanAbsolutePercentageError()`
   * `checkErrorConvergence()`
   * `calculateStatOfMonteCarlo()`
5. Generate visualization outputs:

   * Estimated versus actual signals
   * Convergence and error plots
   * Deviation impact analysis

### Example Code Snippet

```matlab
% In main.m
simulationFileName = 'sim_pf_experiment';
parameter_ekf_ukf_pf; % Load system parameters

% Execute Simulink simulation and analyze results
% Generates plots and .mat files containing SignalAnalysisTool objects
```

## 📦 Project File Structure

```plaintext
.
├── main.m                    # Main simulation and analysis script
├── SignalAnalysisTool.m      # Analysis class definition
├── parameter_ekf_ukf_pf.m    # System and filter parameter configuration
├── proj_sim_*.slx            # Simulink models for EKF, UKF, and PF
├── *.mat                     # Saved simulation results
├── *.txt                     # Simulation log files
```

## 📘 License

This project is licensed under the MIT License. Please refer to the `LICENSE` file for details.

## 👨‍💻 Author

Azka Hariz Sartono
Department of Electrical Engineering
Universitas Gadjah Mada
