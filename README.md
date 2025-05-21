# EVTOL Network Optimization and Demand Estimation

This repository presents a comprehensive framework for the **design and simulation of an EVTOL (Electric Vertical Take-Off and Landing) transportation network**, focusing on demand estimation, infrastructure planning, vehicle modeling, and network optimization.

---

## Project Highlights

### 1. Demand Estimation
- Utilizes **ACS (American Community Survey)** and **LODES (LEHD Origin-Destination Employment Statistics)** data from 2019
- Estimates the **volume and spatial distribution of potential EVTOL demand**
- Spatial processing and aggregation of home-work travel pairs

### 2. Station Location Optimization
- Applies **K-Means Clustering** to identify **optimal station locations**
- Clusters high-demand zones to propose **minimally redundant and spatially balanced hubs**

### 3. EVTOL Modeling – *Spricho*
- Models **Spricho EVTOL** as a **point-mass system**
- Considers aerodynamic constraints, power-to-weight ratio, and energy limits
- Built with a focus on feasibility and system-level performance

### 4. Monte Carlo Network Simulation
- Simulates the **EVTOL transportation network in MATLAB/Simulink**
- Accounts for variability in:
  - Passenger distribution
- Provides **performance metrics** like delay, energy use, and throughput

### 5. Optimization of System Design
- Optimizes the **combination of EVTOL fleet size and station count**
- Objective: **Minimize a cost function** that may include
  - Operational cost
  - Passenger wait time
  - Evacuation Time
  - Waiting Time
  - ...

---

