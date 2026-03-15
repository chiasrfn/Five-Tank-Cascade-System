# Five-Tank-Cascade-System

This repository contains the code, models, and simulation results developed for the **Networked Control** course within the Master's degree program in Engineering.

## Project Overview
The objective of this project is to control a system of **five tanks in cascade** using various state-feedback control structures. The system is modeled in **discrete-time** and evaluated through four different communication and control architectures:

* **Centralized:** A single controller managing all tanks.
* **Decentralized:** Independent controllers for each tank.
* **Distributed (Unidirectional):** A string-based communication scheme where data flows in one direction.
* **Distributed (Bidirectional):** A string-based communication scheme with two-way data exchange.

### Control Techniques
The following control strategies were implemented using Linear Matrix Inequalities (LMIs):
- **Simple LMI:** Basic stabilization.
- **Disk Placement:** LMI with eigenvalue constraints within a specific disk ($\alpha, \rho_{target}$).
- **$H_2$ Control:** Optimization for noise and disturbance rejection.
- **$H_{\infty}$ Control:** Optimization for robust performance against worst-case disturbances.

---

## Repository Structure

### MATLAB Scripts
The core logic is located in the `MATLAB/` directory:

* `IrrigationChannel.m`: **Main Entry Point.** Initializes the system, runs simulations, and allows the selection of different controllers for analysis.
* `LMI_DT_DeDicont.m`: Implements the standard LMI for discrete-time systems.
* `LMI_DT_Disk_Struct.m`: Constructs the LMI for eigenvalue placement within a disk $(\alpha, \rho)$.
* `LMI_DT_H2.m`: Setup for the $H_2$ optimal controller.
* `LMI_DT_Hinf.m`: Setup for the $H_{\infty}$ robust controller.
* `simulate_closed_loop.m`: Handles the closed-loop simulation of the system.
* `di_fixed_modes.m`: Computes the fixed modes of the system to analyze decentralizability.
* `circle.m`: Utility to plot the unit circle for stability and eigenvalue analysis.
### Images
Simulation results are stored in the `images/` directory, organized into two subfolders:

* **`control_action/`**: Contains system response plots. Files are named as `[controller]_[structure]_h[level]`.
* **`eigenvalues/`**: Contains plots of the closed-loop eigenvalues within the unit circle. Files are named as `[controller]_[structure]`.

---

## How to Use
1.  Ensure you have **MATLAB** installed with the **Robust Control Toolbox** or an LMI solver (like SeDuMi or YALMIP).
2.  Clone the repository
3.  Navigate to the `MATLAB/` folder.
4.  Run the main script:
    ```matlab
    IrrigationChannel
    ```
5.  Follow the command window prompts to select the control architecture and the specific controller you wish to analyze.

---

## 🎓 Academic Context
This project was developed as part of the **Networked Control Systems** curriculum. It demonstrates the trade-offs between communication complexity and control performance in interconnected systems.
