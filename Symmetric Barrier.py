import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def process_particle_data(file_path, boundary=848, position_error=15, time_error=1/30):
    # Load the CSV file
    data = pd.read_csv(file_path)

    # Dynamically identify particle columns
    particle_columns = [
        (col, col.replace("X_Particle", "Y_Particle"))
        for col in data.columns if "X_Particle" in col
    ]

    # Initialize lists to store results
    ratios = []
    d_ratios = []
    times = []

    for _, row in data.iterrows():
        N_L = 0  # Particles on the left side
        N_R = 0  # Particles on the right side

        for x_col, y_col in particle_columns:
            x_coord = row[x_col]

            # Determine if the particle is on the left or right side
            if x_coord < boundary:
                N_L += 1
            elif x_coord > boundary:
                N_R += 1

        # Calculate the ratio for the frame
        if (N_R + N_L) > 0:
            ratio = (N_R - N_L) / (N_R + N_L)
            ratios.append(ratio)

            # Error propagation for Ratio
            dN_R = np.sqrt(N_R) * position_error
            dN_L = np.sqrt(N_L) * position_error

            #d_ratio = np.sqrt(
            #    ((1 / (N_R + N_L))**2 * dN_R**2) +
            #    ((1 / (N_R + N_L))**2 * dN_L**2) +
            #    ((N_R - N_L) / (N_R + N_L)**2)**2 * (dN_R**2 + dN_L**2)
            #)

            d_ratio = np.sqrt(
                ((2*N_L*dN_R /(N_R + N_L)**2) ** 2) +
                ((2*N_R * dN_L/ (N_R + N_L)**2) ** 2)

            )
            d_ratios.append(d_ratio)
        else:
            ratios.append(0)
            d_ratios.append(0)

        times.append(row['Time_s'])

    # Convert to  arrays for calculations
    ratios = np.array(ratios)
    d_ratios = np.array(d_ratios)
    times = np.array(times)

    # Calculate mean and variance of ratios
    mean_ratio = np.mean(ratios)
    variance_ratio = np.var(ratios)

    # Error propagation for mean
    mean_error = np.sqrt(np.sum(d_ratios**2)) / len(ratios)


    # Calculate variance error
    #variance_error = np.sqrt(np.sum(2 * ratios * d_ratios)**2) / len(ratios)
    contribution_from_ratios = (1 / len(ratios)) * np.sum(2 * (ratios - mean_ratio) * d_ratios)
    contribution_from_mean = (1 / len(ratios)) * np.sum(2 * (ratios - mean_ratio) * mean_error)

    variance_error = np.sqrt(contribution_from_ratios ** 2 + contribution_from_mean ** 2)

    # Print results
    print(f"Mean Ratio: {mean_ratio} ± {mean_error}")
    print(f"Variance of Ratio: {variance_ratio} ± {variance_error}")

    # Plot the ratio as a function of time
    plt.plot(times, ratios)
    plt.xlabel('t [s]')
    plt.ylabel(r'$\frac{N_R - N_L}{N_R + N_L}$ [unitless]')
    #plt.title('Ratio as a Function of Time')
    plt.grid()
    plt.legend()
    plt.show()

# Choose file 
file_path = r"C:\Users\anna1\PycharmProjects\LAB C - projects\Active Matter 2\sym14.csv"  # Replace with your actual CSV file path
process_particle_data(file_path)

