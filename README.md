# SCIF20002: SEIR model for disease transmission

## Overview:
This project investigates the SEIR model of disease transmission using two approaches:
1. Coupled differential equations
    - The reduced SEIR differential equations are solved as an Initial Value Problem (IVP)
    - Different values of parameters and initival base case are explored
    - The counts of populations in each compartment are plotted over time
2. Stochastic Monte Carlo model
    - A Brownian movement on a lattice simulation
    - Agents (individuals) move randomly locally interacting with one another
    - The solution uses object-oriented programming

## Setup:
This project uses Python.
### Required libraries
- `numpy`
- `matplotlib`
- `random`
- `scipy` for Part 1 (`solve_ivp`)

## File Structure:
The project contains the following main files
### `part1-seir.ipynb`
This notebook contains the deterministic SEIR model from Part 1.

Main contents:
- definition of the reduced SEIR equations
- numerical solution using `scipy.integrate.solve_ivp`
- verification against the baseline case from the assessment brief
- parameter exploration for:
  - different $R_0$ values
  - different $\gamma$ values
  - different $\sigma$ values
  - different initial conditions
- plots of SEIR compartment populations and infected-only comparison plots

### `agent.py`
Defines the `Agent` class used in the Monte Carlo simulation.

Main functionality:
- stores an agent’s lattice coordinates
- stores the agent’s current SEIR state
- updates position when movement occurs
- updates the agent’s state when instructed by the simulation

### `simulation.py`
Defines the `Simulation` class used in the Monte Carlo model.

Main functionality:
- creates the lattice
- initialises agents randomly on empty sites
- checks that the number of agents does not exceed the number of available lattice sites
- finds valid neighbouring positions
- checks whether an agent has infected neighbours
- moves agents stochastically on the lattice
- updates agent states according to:
  - neighbour infection
  - E to I with probability ($\sigma$)
  - I to R with probability ($\gamma$)
- records the compartment counts over time
- saves lattice snapshots at selected Monte Carlo steps
- plots:
  - SEIR population history over Monte Carlo step
  - lattice snapshots showing the spatial distribution of agents

### `run_mcs.ipynb`
This notebook is used to run and explore the Monte Carlo simulation from Part 2.

Main contents:
- baseline simulation matching the assessment brief
- testing and validation checks
- parameter exploration for:
  - different agent densities
  - different $\gamma$ values
  - different $\sigma$ values
- lattice snapshots and population history plots
- qualitative comparison with Part 1

### `README.md`
Project documentation and usage instructions.

### `.gitignore`
Specifies ignored files such as cache folders, temporary notebook files, editor files, and non-submission documents.

## How to run the project
1. ⁠Launch Jupyter Notebook:
```bash 
jupyter notebook part1-seir.ipynb
```
2. Open the notebook in the ⁠ notebooks/ ⁠ folder and run the cells.

## Contribution and git workflow
This project was developed using Git with feature branches.

A typical workflow for collaborators would be:
### 1. Clone the repository
```bash
git clone <repository-url>
cd <repository-folder>
```

### 2. Create a new branch
```bash
git checkout -b your-branch-name
```

### 3. Make changes and check status
```bash
git status
```

### 4. Stage and commit changes
```bash
git add file-name(s)
git commit -m "Describe the change clearly"
```

### 5.Switch branches when needed
```bash
git checkout main
git checkout your-branch-name
```

### 6. Merge completed work into main
```bash
git checkout main
git merge your-branch-name
```

## Possible future extensions
Possible improvements to the current model include:

- allowing reinfection through a small reinfection probability
- modelling weakening immunity over time
- introducing long-distance movement to represent travel
- allowing different transition probabilities for reinfected agents
- adding more realistic movement or contact behaviour