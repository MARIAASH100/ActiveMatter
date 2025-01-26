import pandas as pd
import numpy as np

# Parameters
N_square = 1000  # Number of random points generated per iteration in the square
N_lim = 14 # Target number of particles inside the circle
num_frames = 3600  # Number of frames
circle_center = (537, 971)  # Circle center (x, y)
radius = 435  # Circle radius
frame_time_step = 1 / 30  # Time difference between frames ( 30 FPS)

# Prepare the dataframe
columns = ["Frame", "Time_s", "ParticlesDetected", "X_CenterOfMass", "Y_CenterOfMass"]
data = pd.DataFrame(columns=columns)

# Iterate over frames
for frame_idx in range(num_frames):
    x_points = []
    y_points = []

    # Generate points until we reach N_lim inside the circle
    while len(x_points) < N_lim:
        # Generate random points in the square
        x_square = np.random.uniform(0, 1914, N_square)
        y_square = np.random.uniform(0, 1075, N_square)

        # Filter points inside the circle
        distances = np.sqrt((x_square - circle_center[0])**2 + (y_square - circle_center[1])**2)
        inside_circle = distances <= radius

        # Add valid points to the list
        x_points.extend(x_square[inside_circle])
        y_points.extend(y_square[inside_circle])

    # Trim to exactly N_lim points
    x_points = np.array(x_points[:N_lim])
    y_points = np.array(y_points[:N_lim])

    # Calculate center of mass
    x_com = np.mean(x_points)
    y_com = np.mean(y_points)

    # Prepare frame data
    frame_data = {
        "Frame": frame_idx + 1,
        "Time_s": frame_idx * frame_time_step,
        "ParticlesDetected": N_lim,
        "X_CenterOfMass": x_com,
        "Y_CenterOfMass": y_com,
    }

    # Add particle positions for this frame
    for i in range(N_lim):
        frame_data[f"X_Particle_{i+1}"] = x_points[i]
        frame_data[f"Y_Particle_{i+1}"] = y_points[i]

    # Append to dataframe
    data = pd.concat([data, pd.DataFrame([frame_data])], ignore_index=True)

# Save to CSV:
file_path = "ideal_gas_fixed_circle_particles.csv"
data.to_csv(file_path, index=False)
print(f"CSV saved at {file_path}")
