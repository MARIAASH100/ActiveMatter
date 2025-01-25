# ActiveMatter
We used MATLAB and Python to investigate the motion of active matter modeled by Hexbug robots. The study explores self-alignment in various arenas, the emergence of collective alignment as more Hexbugs are introduced, and the impact of symmetric and asymmetric barriers on their behavior.
## find_pixel_zise.m
This code is applicable for all analyses and is used to define regions or points of interest within a video using pixel-based selection. By manually selecting the top-left and bottom-right corners of the region of interest (ROI), it identifies the area in pixels and allows further calculations for radii or other measurements. Additionally, this process establishes a relationship between pixel dimensions and real-world units, such as centimeters, enabling accurate interpretations and scaling of objects like an arena or other features within the video.
## Hexbug_video_analysis_for_students_single_first.m
This code processes a video to track the trajectory of a single particle (e.g. Hexbug) using the Hough transform for circle detection. To ensure accurate results, the camera should be positioned high enough to capture the Hexbug's movement within its cage while avoiding edge aberrations, which could interfere with detection and tracking. After positioning, the minimum and maximum circle radii (R_min and R_max)  are defined through "find_pixel_size.m" to set the detection bounds for the Hough transform, which identifies dark circular features in each video frame. The program then extracts the particle's X and Y positions frame by frame and saves the trajectory data, including frame numbers, X, and Y coordinates in pixels, into a CSV file named trajectory_data_with_frame_numbers.csv. 







