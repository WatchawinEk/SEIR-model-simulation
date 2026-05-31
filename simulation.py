import numpy as np
import random
from agent import Agent
import matplotlib.pyplot as plt

class Simulation:
    # Constructor
    def __init__(self, box_width, box_height, n_agents, p_susceptible, beta, sigma, gamma):
        # Checks so that inputs are correct
        if n_agents > box_width * box_height:
            raise ValueError(f"Number of agents exceed number of lattice sites {box_width * box_height}")
        
        if box_width <= 0 or box_height <= 0:
            raise ValueError("Grid dimensions must be positive")
        
        if n_agents <= 0:
            raise ValueError("Number of agents must be positive")
        
        if not(0 <= p_susceptible <= 1):
            raise ValueError("Initial susceptible fraction should be between 0 and 1")
        
        if not(0 <= beta <= 1):
            raise ValueError("Beta should be between 0 and 1")
        
        if not(0 <= sigma <= 1):
            raise ValueError("Sigma should be between 0 and 1")
        
        if not(0 <= gamma <= 1):
            raise ValueError("Gamma should be between 0 and 1")

        # Assigning parameters
        self.width = box_width # Width of the grid
        self.height = box_height # Height of the grid
        self.n_agents = n_agents # Total number of agents
        self.p_susceptible = p_susceptible # Fraction of susceptible population at the start
        self.beta = beta # Infection rate
        self.sigma = sigma # Incubation rate
        self.gamma = gamma # Recovery rate

        self.grid = np.zeros((box_height, box_width), dtype=int) # Create grid
        self.agents = [] # To store the agents and their information
        self.initialise_agents()
        self.snapshots = {} # To store snapshot of interested step in MCS

        # To keep track of how many agents are in each state after each step
        self.s_history = []
        self.e_history = []
        self.i_history = []
        self.r_history = []

    # Randomly place the agents on the grid, following the fraction of p_susceptible
    # Careful: Probability of susceptible + Probability of exposed = 1
    def initialise_agents(self):
        visited = 0

        while visited < self.n_agents:
            # Randomly pick out a grid in the box
            x = random.randint(0, self.width - 1)
            y = random.randint(0, self.height - 1)
            
            # Check if it's occupied or not
            # Empty = 0 / S = 1 / E = 2 / I = 3 / R = 4
            if self.grid[y, x] == 0:
                if random.random() < self.p_susceptible:
                    state = 1
                else:
                    state = 2

                agent = Agent(x, y, state) # Initialise an agent instance
                self.agents.append(agent) # Append it to the list of all agents
                self.grid[y, x] = state # Update the grid accordingly
                visited += 1

    # Count the total number of agent in each compartment
    def count_states(self):
        s_count = 0
        e_count = 0
        i_count = 0
        r_count = 0

        for agent in self.agents:
            state = agent.state
            if state == 1:
                s_count += 1
            elif state == 2:
                e_count += 1
            elif state == 3:
                i_count += 1
            elif state == 4:
                r_count += 1
        
        return s_count, e_count, i_count, r_count

    # Count the current stataes and store them as history lists
    def record_history(self):
        s_count, e_count, i_count, r_count = self.count_states()

        self.s_history.append(s_count)
        self.e_history.append(e_count)
        self.i_history.append(i_count)
        self.r_history.append(r_count)

    # To check and verify that, at any step, the total sum of fraction of population 
    # in each compartment is approximate equal to 1 
    def fraction_check(self):
        s_count, e_count, i_count, r_count = self.count_states()

        # The sum of fraction of all compartment should be approximate 1
        frac = (s_count + e_count + i_count + r_count) / self.n_agents 
        
        if abs(1 - frac) <= 0.1:
            return("The sum of fraction is approximately 1")
        else:
            return(f"The sun of fraction {frac} is wrong")

    # Get neighbour's position -> 8 of them
    def get_neighbour(self, x, y):
        neighbours = []

        # Loop through 8 possible neighbour grids
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                
                new_x = x + dx
                new_y = y + dy

                # The wall is hard so only neighbours within the grid is possible
                if 0 <= new_x < self.width and 0 <= new_y < self.height:
                    neighbours.append((new_x, new_y))
        
        return neighbours # A list of neighbours
    
    # Check if the neighbours are infected
    def has_infected_neighbours(self, agent):
        # Loop through the neighbours around and check if they are infected or not
        for x, y in self.get_neighbour(agent.x, agent.y):
            if self.grid[y, x] == 3:
                return True
            
        return False
    
    # Randomly move an agent in the adjacent grids
    # The move only happens if the target cell is within the box and not occupied, otherwise stay put
    def move_agent(self, agent):
        neighbours = self.get_neighbour(agent.x, agent.y) # The check whether it stays within the box happen in get_neighbour
        new_x, new_y = random.choice(neighbours)
        
        # Only move if the grid is not occupied
        if self.grid[new_y, new_x] == 0:
            self.grid[agent.y, agent.x] = 0
            agent.set_position(new_x, new_y) 
            self.grid[new_y, new_x] = agent.state

    # Update the state according to the rate that we inputted
    def update_state(self, agent):
        # Reminder: S = 1 / E = 2 / I = 3 / R = 4
        prev_state = agent.state
        if prev_state == 1:
            if self.has_infected_neighbours(agent): # If at least one neighbour is infecteted, the agent becomes exposed
                agent.set_state(2)
        elif prev_state == 2:
            if random.random() < self.sigma: # E becomes I with probability sigma
                agent.set_state(3)
        elif prev_state == 3:
            if random.random() < self.gamma: # I becomes R with probability gamma
                agent.set_state(4)

        # Also have to update the lattice
        if agent.state != prev_state:
            self.grid[agent.y, agent.x] = agent.state

    # Run 1 Monte Carlo Step (MCS)
    def mcs(self):
        for agent in self.agents:
            self.move_agent(agent)
            self.update_state(agent)
        
        self.record_history()

    # Run the MCS for n_steps
    def run_simulation(self, n_steps, snapshot_steps=None):
        self.record_history() # Initial lattice (without MCS yet)
        self.snapshots = {}

        # If no snapshote steps are provided, use empty list
        if snapshot_steps is None:
            snapshot_steps = []
        else:
            snapshot_steps.sort()

        # Check that all requested steps are within the range we run MCS for
        for snap_step in snapshot_steps:
            if snap_step < 0 or snap_step > n_steps:
                raise ValueError(f"Snapshot step {snap_step} is out of bounds. It must be between 0 and {n_steps}")

        # Save initial lattic if requested
        if 0 in snapshot_steps:
            self.snapshots[0] = self.grid.copy()

        # Run step by step 
        # So i-th step means after i MCS
        for step in range(1, n_steps+1):
            self.mcs()
            
            # Save a copy if the step was requested
            if step in snapshot_steps:
                self.snapshots[step] = self.grid.copy() 
            
            self.record_history()

    # Function to plot the lattice at a chosen step
    def plot_lattice(self, step):
        if step not in self.snapshots:
            raise ValueError(f"No snapshot saved for step {step}.")

        snapshot = self.snapshots[step]

        susceptible_x, susceptible_y = [], []
        exposed_x, exposed_y = [], []
        infected_x, infected_y = [], []
        recovered_x, recovered_y = [], []

        for y in range(snapshot.shape[0]):
            for x in range(snapshot.shape[1]):
                state = snapshot[y, x]

                # Reminder S = 1 / E = 2 / I = 3 / R = 4 
                if state == 1:
                    susceptible_x.append(x)
                    susceptible_y.append(y)
                elif state == 2:
                    exposed_x.append(x)
                    exposed_y.append(y)
                elif state == 3:
                    infected_x.append(x)
                    infected_y.append(y)
                elif state == 4:
                    recovered_x.append(x)
                    recovered_y.append(y)
                else: 
                    continue

        plt.figure(figsize=(7, 7))
        plt.scatter(susceptible_x, susceptible_y, label="Susceptible", s=15)
        plt.scatter(exposed_x, exposed_y, label="Exposed", s=15)
        plt.scatter(infected_x, infected_y, label="Infected", s=15)
        plt.scatter(recovered_x, recovered_y, label="Recovered", s=15)

        plt.xlim(0, self.width)
        plt.ylim(0, self.height)
        plt.xlabel("x position")
        plt.ylabel("y position")
        plt.title("Monte Carlo SEIR simulation")
        plt.legend(loc="upper right")
        plt.show()

    # Function to plot the overall history of each compartment
    def plot_history(self):
        plt.figure(figsize=(9, 5))
        plt.plot(self.s_history, label="Susceptible")
        plt.plot(self.e_history, label="Exposed")
        plt.plot(self.i_history, label="Infected")
        plt.plot(self.r_history, label="Recovered")

        plt.xlabel("Monte Carlo Step")
        plt.ylabel("Population")
        plt.title(f"Monte Carlo SEIR Simulation")
        plt.legend()
        plt.show()