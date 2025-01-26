%% --- Select CSV File ---
%-------------------------
[file, path] = uigetfile({'*.csv', 'CSV Files (*.csv)'}, 'Select the CSV File with Frame, X, Y Data');
if isequal(file, 0)
    disp('No file selected. Exiting...');
    return;
end

% --- Load the CSV File ---
full_path = fullfile(path, file);
data = readmatrix(full_path); % Read the CSV data without headers
frame_numbers = data(:, 1); % Frame numbers (used for time conversion)
x_coords = data(:, 2);      % X-coordinates
y_coords = data(:, 3);      % Y-coordinates

% --- Conversion Parameters ---
pixels_to_cm = 0.06; % Conversion factor: "pixels to cm"=SIZE_cm/SIZE_pixels=>
%Measure the real-world size of the reference object in centimeters SIZE_cm
%Count how many pixels the object spans in the image SIZE_pixels
fps = 30; % Frames per second (adjust based on your data)
frame_interval = 1 / fps; % Time per frame in seconds

% --- Convert Coordinates to Real-World Values ---
%real_distance(cm)=Pixel_distancexpixels_to_cm
x_coords_cm = x_coords * pixels_to_cm; % Convert X-coordinates to cm
y_coords_cm = y_coords * pixels_to_cm; % Convert Y-coordinates to cm

% --- Convert Frame Numbers to Time ---
time_seconds = frame_numbers * frame_interval; % Convert frame numbers to seconds

% --- Calculate MSD ---
N = length(x_coords_cm); % Number of frames
max_lag = min(100, N - 1); % Ensure max_lag does not exceed available frames
msd = zeros(max_lag, 1); % Preallocate MSD array
lag_times = (1:max_lag) * frame_interval; % Lag times in seconds

% Loop over all possible lag times
for lag = 1:max_lag
    displacements = (x_coords_cm(1+lag:end) - x_coords_cm(1:end-lag)).^2 + ...
                    (y_coords_cm(1+lag:end) - y_coords_cm(1:end-lag)).^2;
    msd(lag) = mean(displacements); % Mean squared displacement for this lag
end

%% --- Error Parameters for MSD(lag) ---
%---------------------------------------
Dx = 0.6; % Error in x (cm)
Dy = 0.6; % Error in y (cm)

% Preallocate MSD error array
msd_error = zeros(max_lag, 1);

% Loop over all possible lag times to calculate MSD errors
for lag = 1:max_lag
    x_disp = x_coords_cm(1+lag:end) - x_coords_cm(1:end-lag); % X displacements
    y_disp = y_coords_cm(1+lag:end) - y_coords_cm(1:end-lag); % Y displacements
    
    % Propagate the error for MSD at this lag
    msd_error(lag) = sqrt(sum(4 * (x_disp.^2 * Dx^2 + y_disp.^2 * Dy^2))) / (N - lag);
end

%% --- Error in frame interval <=> time error---
%-----------------------------------------------

%---Changing error---
%--------------------
%frame_interval_error = 1 / fps; % Time error in seconds

% Error in lag times
%lag_time_errors = (1:max_lag)' * frame_interval_error;

%---Constant error---
%--------------------
frame_interval_error = 1 / fps; % Time error in seconds

% Lag time errors (constant for all lags)
lag_time_errors = repmat(frame_interval_error, max_lag, 1); % Same value for all lags

%% --- Plot MSD on Log-Log Scale ---
%-----------------------------------
figure; % Super Important in order for the user to choose points by himself :)
loglog(lag_times, msd, 'o-', 'LineWidth', 1.5);
xlabel('Lag Time (\Delta t) [s]'); % Updated to seconds
ylabel('MSD [cm^2]'); % Updated to cm^2
title('Mean Squared Displacement vs Lag Time (Real-World Units)');
grid on;

disp('MSD calculation with real-world values complete and plot generated.');

% --- Log-Log Transformation ---
%we use it because the upper plot in loglog scaling 
log_lag_times = log10(lag_times); % Log-transform lag times
log_msd = log10(msd); % Log-transform MSD values

% Allow user to select two points => for Ballistic
[x_bal_points, y_bal_points] = ginput(2);

% Sort the selected points to ensure x1<x2
x1_bal  = min(x_bal_points);
x2_bal = max(x_bal_points);

% Define bal segment range based on selected points
bal_idx = log_lag_times >= log10(x1_bal) & log_lag_times <= log10(x2_bal);

% --- Extract 'bal' Segment Points ---
%-----------------------------------
% all the points for the first linear fit(Red)
x_bal = log_lag_times(bal_idx);
y_bal = log_msd(bal_idx);

% Linear Fit for bal Segment
p_bal = polyfit(x_bal, y_bal, 1);
a1_bal = p_bal(1); % Slope
a0_bal = 10^(p_bal(2)); % Intercept (convert back to linear scale)

% Fitted Line for bal Segment
x_fit_bal = linspace(min(x_bal), max(x_bal), 100);
y_fit_bal = p_bal(1) * x_fit_bal + p_bal(2);

% --- Define 'dif' Segment Points ---
%------------------------------------
% all the points for the second linear fit(Blue)
[x_dif_points, y_dif_points] = ginput(2);

% Sort the selected points to ensure x1< x2
x1_dif  = min(x_dif_points);
x2_dif = max(x_dif_points);

% Define dif segment range based on selected points
dif_idx = log_lag_times >= log10(x1_dif) & log_lag_times <= log10(x2_dif);

% --- Extract dif Segment Points ---
x_dif = log_lag_times(dif_idx);
y_dif = log_msd(dif_idx);

% Linear Fit for dif Segment
p_dif = polyfit(x_dif, y_dif, 1);
a1_dif = p_dif(1); % Slope
a0_dif = 10^(p_dif(2)); % Intercept (convert back to linear scale)

% Fitted Line for dif Segment
x_fit_dif = linspace(min(x_dif), max(x_dif), 100);
y_fit_dif = p_dif(1) * x_fit_dif + p_dif(2);

% --- Final Plot with Fitted Functions ---
%-----------------------------------------
% This section will give the plot of MSD+Bal+Dif
figure;
loglog(lag_times, msd, 'o-', 'LineWidth', 1.5); % Original MSD vs Lag Times
hold on;

% Plot fitted functions on the same graph
fitted_bal = 10.^(a1_bal * log10(lag_times) + log10(a0_bal));
fitted_dif = 10.^(a1_dif * log10(lag_times) + log10(a0_dif));

loglog(lag_times , fitted_bal, 'r-', 'LineWidth', 1.5); % bal segment fit
loglog(lag_times, fitted_dif, 'b-', 'LineWidth', 1.5); % dif segment fit

xlabel('\Delta t [s]');
ylabel('MSD [cm^2]');
legend('Data', 'Fit: Short Time Scales (Ballistic Motion)', 'Fit: Long Time Scales (Diffusive Motion)');
grid on;
hold off;

%% --- Constants + First step for error caculation ---
%-----------------------------------------------------

ln10 = log(10); % logarithm of 10

% Logarithmic Transformation and Indirect Errors
log_lag_times = log10(lag_times); % <3
lag_time_errors = lag_time_errors(:); % Ensure column vector=> believe us it it important
lag_times = lag_times(:); % Ensure column vector
log_lag_times_error = lag_time_errors ./ (lag_times * ln10); % Propagate error for log10(lag_times) <3

log_msd = log10(msd); % <3
log_msd_error = msd_error ./ (msd * ln10); % Propagate error for log10(MSD) <3

% Ensure all the relevant points denoted as <3 for the fir is column vector
log_lag_times = log_lag_times(:);  
log_lag_times_error = log_lag_times_error(:);  
log_msd = log_msd(:);  
log_msd_error = log_msd_error(:); 

%% ---Weighted Fits and Error Calculations + the fit parameters---
%-----------------------------------------------------------------

% --- bal(ballistic) Segment Weighted Fit ---
%--------------------------------------------
x_bal = log_lag_times(bal_idx);
y_bal = log_msd(bal_idx);
y_err_bal = log_msd_error(bal_idx); 
x_err_bal = log_lag_times_error(bal_idx); 

% Weights for bal
weights_bal = 1 ./ y_err_bal.^2;

% Fit using polyfit with weights for bal
[p_bal, S_bal] = polyfit(x_bal, y_bal, 1);

a1_bal = p_bal(1); % Slope
a0_bal = p_bal(2); % Intercept

% Error Estimates for bal
cov_matrix_bal = inv(S_bal.R) * inv(S_bal.R)'; % Covariance matrix
Da1_bal = sqrt(cov_matrix_bal(1, 1)); % Uncertainty in slope
Da0_bal = sqrt(cov_matrix_bal(2, 2)); % Uncertainty in intercept

% --- dif(diffusive) Segment Weighted Fit ---
%--------------------------------------------
x_dif = log_lag_times(dif_idx);
y_dif = log_msd(dif_idx);
y_err_dif = log_msd_error(dif_idx);  
x_err_dif = log_lag_times_error(dif_idx);  

% Weights for dif
weights_dif = 1 ./ y_err_dif.^2;

% Fit using polyfit with weights for dif
[p_dif, S_dif] = polyfit(x_dif, y_dif, 1);

a1_dif = p_dif(1); % Slope
a0_dif = p_dif(2); % Intercept

% Error Estimates for dif
cov_matrix_dif = inv(S_dif.R) * inv(S_dif.R)'; % Covariance matrix
Da1_dif = sqrt(cov_matrix_dif(1, 1)); % Uncertainty in slope
Da0_dif = sqrt(cov_matrix_dif(2, 2)); % Uncertainty in intercept

%% --- Display Results---
%------------------------
disp('bal Segment Fit Results:');
fprintf('Slope (a1) = %.4f ± %.4f\n', a1_bal, Da1_bal);
fprintf('Intercept (a0) = %.4f ± %.4f\n', a0_bal, Da0_bal);

disp('dif Segment Fit Results:');
fprintf('Slope (a1) = %.4f ± %.4f\n', a1_dif, Da1_dif);
fprintf('Intercept (a0) = %.4f ± %.4f\n', a0_dif, Da0_dif);

%% ---Intersection Point---
%--------------------------
log_t = (a1_bal - a1_dif) / (a0_dif - a0_bal);

% Error in intersection point
d_log_t = sqrt( ...
    (1 / (a1_dif - a1_bal))^2 * Da0_bal^2 + ...
    (-1 / (a1_dif - a1_bal))^2 * Da0_dif^2 + ...
    ((a0_bal - a0_dif) / (a1_dif - a1_bal)^2)^2 * Da1_bal^2 + ...
    (-(a0_bal - a0_dif) / (a1_dif - a1_bal)^2)^2 * Da1_dif^2 ...
);

t = 10^log_t;
dt = abs(d_log_t/(ln10 * log_t));

D_theta = 1 / t;
d_D_theta = dt / t^2;

disp('Intersection Point Results:');
fprintf('t = %.4f ± %.4f\n', t, dt);
fprintf('D_theta = %.4f ± %.4f\n', D_theta, d_D_theta);

