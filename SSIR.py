import numpy as np
import random
import matplotlib.pyplot as plt


def main():
    R0 = 1.34
    tau_i = 6.3
    population_size = 1000
    t_end = 200
    N = 100 # N = number of simulations

    
    final_S = []
    final_I = []
    final_R = []
    final_t = []  # Array to store final times
    
    extinction_counter = 0.0 # Counter of Stochastic Extinction Events

    for n in range(N):
        simulation = SIR(R0, tau_i, population_size, t_end)
        simulation.SimulationLoop()
        S, I, R, t = simulation.final_state()

        final_S.append(S)
        final_I.append(I)
        final_R.append(R)
        final_t.append(t)

        if S >= 500:
            extinction_counter += 1

    probability_of_extinction = (extinction_counter / N) * 100

    # Print final data for all N simulations
    print("Final data after each simulation:")
    
    for i in range(N):
        print(f"Simulation {i+1}: S = {final_S[i]}, I = {final_I[i]}, R = {final_R[i]}, Time = {final_t[i]:.2f}")
    
    # Print the probability of stochastic extinction
    print(f"The Extinction Probability is:  {probability_of_extinction}%")
    
    # Create the time evolution plot
    plt.figure(figsize=(8, 6))
    for _ in range(N):
        simulation = SIR(R0, tau_i, population_size, t_end)
        simulation.SimulationLoop()
        simulation.plot_results()
    
    plt.xlabel('Time (Days)')
    plt.ylabel('Population')
    plt.legend(['Susceptible', 'Infected', 'Recovered'], title="Population")
    plt.title(f'Time Evolution of Disease Spread Over {N} Simulations')
    plt.show()

    # Create the final state scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(range(N), final_S, color='blue', label='Final Susceptible')
    plt.scatter(range(N), final_I, color='red', label='Final Infected')
    plt.scatter(range(N), final_R, color='green', label='Final Recovered')
    plt.scatter(range(N), final_t, color='purple', label='Final Time', marker='x')

    plt.xlabel('Simulation Run')
    plt.ylabel('Population Count')
    plt.legend()
    plt.title(f'Final State of S, I, R and Time After {N} Simulations')
    plt.show()


class SIR:

    def __init__(self, R0, tau_i, population_size, t_end):
        # Parameters
        self.R0 = R0
        self.tau_i = tau_i
        self.beta = R0 / tau_i
        self.gamma = 1 / tau_i
        self.scaled_beta = self.beta * .01

        # Settings
        self.population_size = population_size
        self.t_end = t_end

        # Populations
        self.S = [population_size - 1]
        self.I = [1]
        self.R = [0]
        self.t = [0]

    def propensities(self, S, I):
        # propensities calc
        N_pop = S + I + self.R[-1]
        effective_beta = self.scaled_beta * S / N_pop
        Propensity_I = effective_beta * S * I
        Propensity_R = self.gamma * I
        return Propensity_I, Propensity_R

    def population_update(self, S, I, R, Propensity_I, Propensity_R):
        Propensities = [Propensity_I, Propensity_R]
        Propensity_Sum = Propensities[0] + Propensities[1]

        if Propensity_Sum == 0:
            return S, I, R, 0
        
        # time step
        time_step = np.random.exponential(scale=1 / Propensity_Sum)
        rand_num = random.uniform(0, 1)

        # Event Picker
        if rand_num * Propensity_Sum <= Propensity_I and S > 0:
            new_infected = min(S, 1)  # Control max infections per step
            S -= new_infected
            I += new_infected
        elif I > 0:
            I -= 1
            R += 1

        return S, I, R, time_step
    
    def SimulationLoop(self):
        while self.t[-1] < self.t_end and self.S[-1] + self.I[-1] + self.R[-1] > 0:
            S, I, R = self.S[-1], self.I[-1], self.R[-1]

            if S == 0 or I == 0:
                break

            # Get the propensities and update populations
            propensity_I, propensity_R = self.propensities(S, I)
            S, I, R, time_step = self.population_update(S, I, R, propensity_I, propensity_R)

            if time_step == 0:
                break
            
            # Update populations and time
            self.S.append(S)
            self.I.append(I)
            self.R.append(R)
            self.t.append(self.t[-1] + time_step)

    def plot_results(self):
        # Plot the evolution of the populations over time for this simulation
        plt.plot(self.t, self.S, color='blue', alpha=0.5)
        plt.plot(self.t, self.I, color='red', alpha=0.5)
        plt.plot(self.t, self.R, color='green', alpha=0.5)

    def final_state(self):
        return self.S[-1], self.I[-1], self.R[-1], self.t[-1]
    
if __name__ == "__main__":
    main()