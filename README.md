# ActiveMatter
We used MATLAB and Python to investigate the motion of active matter modeled by Hexbug robots. The study explores self-alignment in various arenas, the emergence of collective alignment as more Hexbugs are introduced, and the impact of symmetric and asymmetric barriers on their behavior.
## Errors
The error estimation in our analysis accounts for Gaussian errors in the determination of particle center points, as particles are identified by fitting circles around them. The radius of these circles serves as the basis for error estimation, and we chose one standard deviation of the radius as the error in both Δx and Δy =>  Δx = Δy = 1σ f the radius. While this approach might slightly overestimate the errors, further refinement and validation are needed. All additional errors were handled using standard indirect error propagation methods to ensure consistency in the analysis. The error in time was estimated based on the frame rate of the video(fps) Δt=1/fps. 
fps), with the time error taken as the reciprocal of the frame rate 
## find_pixel_zise.m
This code is applicable for all analyses and is used to define regions or points of interest within a video using pixel-based selection. By manually selecting the top-left and bottom-right corners of the region of interest (ROI), it identifies the area in pixels and allows further calculations for radii or other measurements. Additionally, this process establishes a relationship between pixel dimensions and real-world units, such as centimeters, enabling accurate interpretations and scaling of objects like an arena or other features within the video.
In this context, the columns represent the X-coordinates, and the rows represent the Y-coordinates, with the region starting from the top-left corner of the frame, designated as the origin (0, 0).
## Hexbug_video_analysis_for_students_single_first.m (by D.Shohat)
This code processes a video to track the trajectory of a single particle (e.g. Hexbug) using the Hough transform for circle detection. To ensure accurate results, the camera should be positioned high enough to capture the Hexbug's movement within its cage while avoiding edge aberrations, which could interfere with detection and tracking. After positioning, the minimum and maximum circle radii (R_min and R_max)  are defined through "find_pixel_size.m" to set the detection bounds for the Hough transform, which identifies dark circular features in each video frame. The program then extracts the particle's X and Y positions frame by frame and saves the trajectory data, including frame numbers, X, and Y coordinates in pixels, into a CSV file named trajectory_data_with_frame_numbers.csv. 
## Hexbug_video_analysis_for_students_green_stickers_first.m (by D.Shohat)
This code tracks green stickers in a video using the Hough transform to detect circular features. To ensure accurate tracking, the camera should be positioned lower to closely capture the arena. Key parameters to configure include N_particles(number of particles in the video), R_min and R_max (minimum and maximum circle radii in pixels, determined from the "find_pixel_zise.m"), Pixel_to_cm (conversion factor based on a known size, like the arena or sticker diameter), and fps (frames per second of the video).
The output includes data columns: Frame	Time_s	X_Particle_1	Y_Particle_1	X_Particle_2	Y_Particle_2...X_Particle_(N_particles)	Y_Particle_(N_particles)	ParticlesDetected	X_CenterOfMass	Y_CenterOfMass.
All tracking data is saved in a CSV file named particle_tracking_data.csv in the same directory as the input video.
Note: This green sticker tracking script uses the frame_isolate_green_try2 function to isolate green regions in each video frame by applying color thresholding in the HSV color space, ensuring accurate detection of the stickers.
## Heatmap_Position.m
This code generates a heatmap of particle positions from a CSV file by aggregating their occurrences across frames, which was used in parts 1 and 3 to visualize self-alignment, flocking, and particle gathering under specific conditions.
# Part 1 - Single Particle
# MSD_lag_single.m
## Oriantential_Coralleation.m
This code calculates the orientational correlation function C(Δt) over time lags for particles, including error propagation, identifies peaks in C(Δt) and fits an exponential decay model to estimate the rotational diffusion constant D_θ with associated uncertainties.

# Part 2 - Circular Arena - Flocking in a Group of Hexbugs
In this section, we explore the phenomenon of flocking behavior among active particles, examining how their collective dynamics deviate from randomness. Additionally, we compare these observations with an ideal gas simulation to highlight the unique features of coordinated motion in active systems. This was achieved using the pair correlation function, an excellent statistical tool for analyzing spatial organization and detecting patterns of clustering or uniformity in particle distributions.
## Ideal Gas Monte Carlo sim.py
oundary conditions play a crucial role in determining particle behavior. In confined systems, the distribution of particles, often analyzed via the pair correlation function g(r), can vary significantly depending on the geometry. For example, the behavior of particles in a circular or spherical confinement (as in an electron gas or Fermi gas) often differs from that in a rectangular box. Simulations and experiments have shown that spherical boundaries can lead to deviations in g(r), reflecting how confinement geometry influences particle interactions and spatial distributions. Because of that Monte Carlo simulation implemented here models a system of particles confined within a circular area to mimic experimental conditions. Specifically, the number of particles (N_lim) corresponds to the number of hexbugs observed in the experiment, ensuring equivalent conditions. The simulation also replicates the same spatial dimensions as the experimental setup to maintain consistency. A total of 3600 frames is generated, matching the duration and frame count of the hexbug video, which is essential for achieving robust statistical analysis. The simulation works by randomly generating points within a larger rectangular area and filtering those that fall within the circular region of interest. This process continues until the desired number of particles is reached for each frame. The center of mass of the particles is then calculated and recorded, along with their individual positions, for each frame.
## Pair Corralation generate CSV file general.py
The calculation of g(r) began with extracting particle positions from each frame of the simulation or experimental data. All pairwise distances between particles were calculated, excluding self-pairs. These distances were binned into a histogram with a bin width of 
Δr=30 pixels, representing the minimum realistic distance between particles. The histogram was normalized by the total number of particle pairs to create the probability density function (PDF) of distances. This PDF was then normalized by the shell area for each bin to account for the geometry of the system, converting it to the local pair density. Finally, the local pair density was normalized by the system-wide average density to produce the pair correlation function g(r). his process was repeated for all frames, and the final g(r) was obtained by averaging the results over all frames to ensure statistical reliability. This approach capturing the physical behavior of both the ideal gas and Hexbugs systems under consistent conditions.
## Pair Correlation Function Hexbugs vs Ideal Gas.py
This script plots the pair correlation functions (PCFs) of Hexbugs and the Ideal Gas on the same graph for direct comparison.

Note: Using consistent file naming conventions like gr_vs_r_14p for Hexbugs and gr_vs_r_ideal_gas14p for the Ideal Gas, with the particle count indicated in the filename => makes life easier and happier ! 

# Part 3 - Group of Hexbugs and how Symmetric and Asymmetric Barrier Impact Them
The reason for studying Hexbugs with symmetric and asymmetric barriers is to observe how particle alignment and collective behavior influence their movement. In a symmetric barrier, the goal is to see if one particle crossing the barrier triggers self-alignment, encouraging others to follow and pass through as a group. However, in an asymmetric barrier, the asymmetry introduces a directional bias, making it more likely for particles to favor one side (typically the left in our setup). This bias disrupts the natural alignment observed in symmetric setups, leading to an uneven particle distribution and unique dynamics across the barrier.
## Symmetric Barrier.py
This code analyzes particle data to calculate the ratio of particles on the left and right sides of a defined barrier over time. The barrier's boundary is set using a "find_pixel_zise.m" (a constant line at "boundary"), enabling the code to count the particles on each side and compute the time-dependent ratio and its associated errors.
## A-Symmetric Barrier.py
This code analyzes particle data in an asymmetric barrier configuration where the asymmetry causes more particles to accumulate on the left side, leading to a mean ratio of less than 0. The particle distribution is determined using three distinct regions: two defined by linear boundaries (derived using "find_pixel_zise.m") and one by a constant vertical line. These boundaries divide the area into left and right regions, enabling the calculation of particle counts and the time-dependent ratio, , along with its associated statistical errors.

# Part 4 - Movable Obstacles
## moving_obstacles.m
This code analyzes the movement of obstacles in the presence of Hexbugs by tracking their displacement relative to an initial reference point over time. The user begins by selecting a video file and marking an initial reference point in the first frame, such as the center of mass of the obstacle. The code then processes every frameJump frame, where we can manually select the new position of interest in each frame. The Euclidean distance between the selected point and the reference point is calculated and converted to cm using a predefined scaling factor. These distances, along with the corresponding time stamps, are stored and used to plot the displacement over time. For symmetrical obstacles the center of mass is expected to remain relatively stable, indicating minimal movement due to the lack of directional bias. In contrast, asymmetrical obstacles cause Hexbugs to gather in some areas, resulting in a noticeable shift in the center of mass over time. This shift is reflected in the distance vs time plot and highlights how obstacle asymmetry influences the collective behavior of the Hexbugs.












