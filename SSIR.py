import numpy as np
import random
import matplotlib.pyplot as plt

# Parameters for the simulation
R0 = 2.5  # Basic Reproduction Number 
tau_i = 10  # Infectious Period (days)
beta = R0 / tau_i  # Transmission rate
gamma = 1 / tau_i  # Recovery rate

scaled_beta = beta * 0.01  # Scaling factor to adjust infection speed

# Number of simulations to run
N = 10  

# Colors for S, I, R populations
color_S = 'blue'
color_I = 'red'
color_R = 'green'

# Run the simulation N times
for n in range(N):
    # Initialize population values
    S = [999]  # Initial susceptible population
    I = [1]    # Initial infected individuals
    R = [0]    # Initial recovered individuals
    t = [0]    # Time starts at 0
    t_end = 200  # t in days

    while t[-1] < t_end and S[-1] + I[-1] + R[-1] > 0:  # Ensure populations don't go negative
        N_pop = S[-1] + I[-1] + R[-1]

        # Check if there are still susceptible or infected individuals
        if S[-1] == 0 or I[-1] == 0:
            break

        # Adjust transmission rate based on current susceptible population
        effective_beta = scaled_beta * S[-1] / N_pop  # Adjust transmission rate with current susceptible population

        # Propensities
        Propensity_I = effective_beta * S[-1] * I[-1]
        Propensity_R = gamma * I[-1]

        if Propensity_I + Propensity_R == 0:  # Prevent division by zero
            break

        # Propensity array
        Propensities = [Propensity_I, Propensity_R]
        Propensity_Sum = Propensities[0] + Propensities[1]

        # Smaller time steps
        time_step = np.random.exponential(scale=1/Propensity_Sum)
        rand_num = random.uniform(0, 1)

        # Infection step, but limit how many infections can occur
        if rand_num * Propensity_Sum <= Propensity_I and S[-1] > 0:
            # Limit the number of new infections based on available susceptible individuals
            new_infected = min(S[-1], 1)  # Adjust this "1" to control max infections per step
            S.append(S[-1] - new_infected)
            I.append(I[-1] + new_infected)
            R.append(R[-1])
        elif I[-1] > 0:  # Recovery happens
            S.append(S[-1])
            I.append(I[-1] - 1)
            R.append(R[-1] + 1)

        t.append(t[-1] + time_step)

    # Plot the populations for all simulations with color-coding for S, I, and R
    plt.plot(t, S, color=color_S)
    plt.plot(t, I, color=color_I)
    plt.plot(t, R, color=color_R)

# Label the axes and show the plot
plt.xlabel('Time (Days)')
plt.ylabel('Population')
plt.legend(['Susceptible', 'Infected', 'Recovered'], title="Population")
plt.title('Multiple Simulations of Disease Spread (S, I, R)')
plt.show()