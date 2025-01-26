import pandas as pd
import matplotlib.pyplot as plt

# File paths for the two datasets of CSV
file_hexbugs = r"C:\Users\anna1\PycharmProjects\LAB C - projects\Active Matter 2\gr_vs_r_14p.csv" # Replace with your Hexbugs g(r) CSV file path
file_ideal_gas = r"C:\Users\anna1\PycharmProjects\LAB C - projects\Active Matter 2\gr_vs_r_ideal_gas14p.csv" # Replace with your Ideal Gas g(r) CSV file path

# Load the data
data_hexbugs = pd.read_csv(file_hexbugs)
data_ideal_gas = pd.read_csv(file_ideal_gas)

# Extract r and g(r) values
r_hexbugs = data_hexbugs["r (distance in pixels)"]
gr_hexbugs = data_hexbugs["g(r)"]

r_ideal_gas = data_ideal_gas["r (distance in pixels)"]
gr_ideal_gas = data_ideal_gas["g(r)"]

# Plot the data
#plt.figure(figsize=(8, 6))
#plt.plot(r_hexbugs, gr_hexbugs, label="Hexbugs", linestyle='-', marker='o', color='blue')
#plt.plot(r_ideal_gas, gr_ideal_gas, label="Ideal Gas", linestyle='--', marker='s', color='green')

#plt.xlabel("r (distance in pixels)")
#plt.ylabel("g(r)")
#plt.title("Comparison of Pair Correlation Functions ideal gas and Hrxbugs")
#plt.legend()
#plt.grid(alpha=0.5)

# Show the plot
plt.show()

# Create a step-like plot to mimic a histogram outline
    #plt.step(r_values, gr_avg, where='mid', color='blue', linewidth=1.5, label="g(r)")
plt.step(r_hexbugs, gr_hexbugs, label="Hexbugs", where='mid', color='blue', linewidth=1.5)
plt.step(r_ideal_gas, gr_ideal_gas, label="Ideal Gas", where='mid', color='green', linewidth=1.5)


plt.xlabel("r [cm]")
plt.ylabel("g(r) [unitless]")
plt.title("N = 14 particles")
plt.legend()
plt.grid(alpha=0.5)

plt.show()

plt.step(r_hexbugs, gr_hexbugs, label="Hexbugs", where='mid', color='blue', linewidth=1.5)
plt.step(r_ideal_gas, gr_ideal_gas, label="Ideal Gas", where='mid', color='green', linewidth=1.5)


plt.xlabel("r [cm]")
plt.ylabel("g(r) [unitless]")
plt.legend()
plt.grid(alpha=0.5)

plt.show()
