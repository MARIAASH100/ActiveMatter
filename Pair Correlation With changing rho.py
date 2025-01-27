import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculate_gr(file_path, Dr, max_r, arena_area, fps, pixel_to_cm):
    # Load the CSV file
    data = pd.read_csv(file_path)

    # Dynamically identify particle columns
    particle_columns = [
        (col, col.replace("X_Particle", "Y_Particle"))
        for col in data.columns if "X_Particle" in col
    ]

    # Define parameters
    r_bins = np.arange(Dr, max_r + Dr, Dr)  # Bins: Dr, 2Dr, ..., max_r
    gr_frames = []  # Store g(r) for each frame

    # Iterate over each frame
    for frame_idx in range(len(data)):
        # Step 1: Extract particle positions for the frame
        positions = []
        for x_col, y_col in particle_columns:
            x = data.loc[frame_idx, x_col]
            y = data.loc[frame_idx, y_col]
            if x != 0 and y != 0:  # Only include detected particles
                positions.append([x, y])
        positions = np.array(positions)
        num_particles = len(positions)

        if num_particles < 2:
            continue  # Skip frames with less than 2 particles

        # Calculate rho_0 dynamically based on ParticlesDetected
        num_particles_real = data.loc[frame_idx, "ParticlesDetected"]  # Get ParticlesDetected from CSV
        rho_0 = num_particles_real / arena_area  # Update rho_0 for each frame

        # Step 2: Calculate all pairwise distances d_ij
        distances = []
        for i in range(num_particles):
            for j in range(num_particles):
                if i != j:  # Avoid self-pairs
                    d_ij = np.sqrt((positions[i][0] - positions[j][0])**2 + (positions[i][1] - positions[j][1])**2)
                    distances.append(d_ij)
        distances = np.array(distances)

        # Step 3: Bin distances into histogram
        hist, _ = np.histogram(distances, bins=r_bins)

        # Step 4: Normalize by the total number of pairs
        total_pairs = num_particles * (num_particles - 1)  # N * (N-1) since distances include both directions
        pdf_dij = hist / total_pairs  # PDF of d_ij

        # Step 5: Normalize by shell area
        shell_areas = 2 * np.pi * r_bins[:-1] * Dr  # A(r) = 2*pi*r*Dr
        num_real_pairs = num_particles * (num_particles - 1) / 2
        normalized_pdf = num_real_pairs * pdf_dij / shell_areas  # Local pair density

        # Step 6: Normalize by system-wide average density
        gr = normalized_pdf / rho_0  # Pair correlation function

        # Store g(r) for this frame
        gr_frames.append(gr)

    # Step 7: Average g(r) over all frames
    gr_avg = np.mean(gr_frames, axis=0)

    # Step 8: Plot g(r) as a function of r
    r_values = r_bins[:-1] + Dr / 2  # Bin centers
    r_values_cm = r_values * pixel_to_cm  # Convert r-values to centimeters

    plt.figure(figsize=(8, 6))
    plt.plot(r_values_cm, gr_avg, label="g(r)", marker='o')
    plt.xlabel("r [cm]")
    plt.ylabel("g(r)")
    plt.title("Pair Correlation Function")
    plt.legend()
    plt.grid()
    plt.show()

    # Save g(r) and r_values to a CSV file
    gr_data = pd.DataFrame({
        "r (distance in cm)": r_values_cm,
        "g(r)": gr_avg
    })
    output_path = "gr_vs_r.csv"
    gr_data.to_csv(output_path, index=False)
    print(f"g(r) data saved to {output_path}")

    return r_values_cm, gr_avg

# Example usage
file_path = r"C:\\Users\\anna1\\PyCharmProjects\\LAB C - projects\\Active Matter 2\\13p.csv"
Dr = 30  # Bin width
max_r = 900  # Maximum distance
arena_area = np.pi * 435**2  # Arena area
fps = 1 / 30  # Frame sampling period in seconds
pixel_to_cm = 1 / 29  # Conversion factor: 1 pixel = 1/29 cm

r_values, gr_avg = calculate_gr(file_path, Dr, max_r, arena_area, fps, pixel_to_cm)
