import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def linear_line_x_of_y(x1, y1, x2, y2):
    # Calculate slope
    m = (x2 - x1) / (y2 - y1)
    # Calculate intercept
    b = x1 - m * y1
    # Return the equation as a function
    def x_of_y(y):
        return m * y + b
    return x_of_y

def process_particle_data_asymmetric(file_path, position_error=15):
    # Define the points for the asymmetric barrier (use pixel_fins.m script)
    r1 = (1074, 148)
    r2 = (805, 444)
    r3 = (805, 619)
    r4 = (1074, 895)

    # Define linear functions for the boundaries
    X_linear1 = linear_line_x_of_y(r1[0], r1[1], r2[0], r2[1])
    X_linear2 = linear_line_x_of_y(r3[0], r3[1], r4[0], r4[1])

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
            y_coord = row[y_col]

            # Check region and count particles
            if 0 < y_coord <= r2[1]:
                if x_coord < X_linear1(y_coord):
                    N_L += 1
                else:
                    N_R += 1
            elif r2[1] < y_coord <= r3[1]:
                if x_coord < r2[0]:
                    N_L += 1
                else:
                    N_R += 1
            elif r3[1] < y_coord <= r4[1]:
                if x_coord < X_linear2(y_coord):
                    N_L += 1
                else:
                    N_R += 1

        # Calculate the ratio for the frame
        if (N_R + N_L) > 0:
            ratio = (N_R - N_L) / (N_R + N_L)
            ratios.append(ratio)

            # Error propagation for Ratio
            dN_R = np.sqrt(N_R) * position_error
            dN_L = np.sqrt(N_L) * position_error

            d_ratio = np.sqrt(
                ((1 / (N_R + N_L))**2 * dN_R**2) +
                ((1 / (N_R + N_L))**2 * dN_L**2) +
                ((N_R - N_L) / (N_R + N_L)**2)**2 * (dN_R**2 + dN_L**2)
            )
            d_ratios.append(d_ratio)
        else:
            ratios.append(0)
            d_ratios.append(0)

        times.append(row['Time_s'])

    # Convert to numpy arrays for calculations
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
    #plt.title('Ratio as a Function of Time (Asymmetric Barrier)')
    plt.grid()
    plt.legend()
    plt.show()

# Choose file
file_path = r"C:\Users\anna1\PycharmProjects\LAB C - projects\Active Matter 2\midAS14second.csv"  # Replace with your actual CSV file name
process_particle_data_asymmetric(file_path)
