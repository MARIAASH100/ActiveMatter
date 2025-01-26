%% ---Select and Read CSV File---
%--------------------------------
[file, path] = uigetfile('*.csv', 'Select a CSV File');
if isequal(file, 0)
    disp('No file selected. Exiting...');
    return;
end

% --- Load the CSV File ---
full_path = fullfile(path, file);
data = readmatrix(full_path); % Read the CSV data without headers
frames = data(:, 1); % Frame numbers (used for time conversion)
X = data(:, 2);      % X-coordinates
Y = data(:, 3);      % Y-coordinates

%% Initialize Variables
num_frames = length(frames);
theta = zeros(num_frames - 1, 1); % Orientation angles in radians
n_vectors = zeros(num_frames - 1, 2); % Orientation vectors [cos(theta), sin(theta)]

%% Calculate Orientations and Orientation Vectors with errors

% Initialize errors
DX = 0.6; % Error in X-coordinate
DY = 0.6; % Error in Y-coordinate

theta_error = zeros(num_frames - 1, 1); % Preallocate theta error
cos_error = zeros(num_frames - 1, 1); % Preallocate cos(theta) error
sin_error = zeros(num_frames - 1, 1); % Preallocate sin(theta) error

for j = 1:num_frames - 1
    dx = X(j + 1) - X(j);
    dy = Y(j + 1) - Y(j);
    r_squared = dx^2 + dy^2; % Distance squared

    % Errors in dx and dy
    dx_error = sqrt(2) * DX;
    dy_error = sqrt(2) * DY;

    theta(j) = atan2(dy, dx); % Compute angle using atan2
    % Error in theta
    theta_error(j) = sqrt(((dy / r_squared)^2 * dx_error^2) + ((dx / r_squared)^2 * dy_error^2));

    n_vectors(j, :) = [cos(theta(j)), sin(theta(j))]; % Compute orientation vector

    % Errors in orientation vector components
    cos_error(j) = abs(sin(theta(j))) * theta_error(j); % Error in cos(theta)
    sin_error(j) = abs(cos(theta(j))) * theta_error(j); % Error in sin(theta)
end

% Combine errors into the same format as n_vectors
n_vector_errors = [cos_error, sin_error]; % Combine errors into a 2D array


 
%% Compute Correlation Function
max_lag = 100; % Maximum frame lag to compute
C = zeros(max_lag, 1); % Correlation values
C_error = zeros(max_lag, 1); % Correlation errors
fps = 30; % Frames per second (adjust if needed)
delta_t = 1 / fps; % Time interval between frames

for n = 1:max_lag
    correlation_sum = 0; % Sum of correlations for the current lag
    error_sum = 0; % Sum of squared errors for the current lag =>for each iteration a new error sum implimented sum(error of n_t dot n_t+Dt)^2
    num_pairs = 0; % Number of valid pairs

    for t = 1:(num_frames - 1 - n)
        % Compute dot product of orientation vectors separated by lag n
        correlation_sum = correlation_sum + dot(n_vectors(t, :), n_vectors(t + n, :));
        % Propagate error for dot product
        dot_product_error = sqrt( ...
            (n_vectors(t + n, 1) * n_vector_errors(t, 1))^2 + ...
            (n_vectors(t, 1) * n_vector_errors(t + n, 1))^2 + ...
            (n_vectors(t + n, 2) * n_vector_errors(t, 2))^2 + ...
            (n_vectors(t, 2) * n_vector_errors(t + n, 2))^2 ...
        );
        % Accumulate squared errors
        error_sum = error_sum + dot_product_error^2; %sum(error of n_t dot n_t+Dt)^2
        num_pairs = num_pairs + 1;
    end

    % Average correlation for the current lag
    if num_pairs > 0
        C(n) = correlation_sum / num_pairs;
        C_error(n) = sqrt(error_sum) / num_pairs;
    end
end
lag_times = (1:max_lag) * delta_t; % Time lags in seconds


%% Error in Lag Times
fps = 30; % Frames per second
delta_t = 1 / fps; % Time per frame
err_t = 1 / fps; % Error in time per frame

% Lag times and their errors
lag_times = (1:max_lag)' * delta_t; % Lag times in seconds
lag_time_errors = (1:max_lag)' * err_t; % Error in lag times


%% ---Identify Peaks---
%----------------------
[pks, locs] = findpeaks(C, lag_times, 'MinPeakProminence', 0.1); % Find peaks 

% Perform logarithmic transformation on the peaks => make like easier with
% the fit
log_pks = log(pks); %easier to impliment polypit
% Convert locs (lag times) to their corresponding indices
[~, locs_indices] = ismember(locs, lag_times);

% Assign errors to the peaks
pks_error = C_error(locs_indices); % Errors in C corresponding to the peaks
log_pks_error = pks_error./log_pks; % Errors in C corresponding to the peaks
lag_times_error_peaks = lag_time_errors(locs_indices); % Errors in lag_times corresponding to the peaks


%% ---Fit Exponential Decay to Peaks using log_pks_error and log_lag_time_error ---
%----------------------------------------------------------------------------------

% Perform weighted linear fit using locs and log_pks
%[fit_params, S] = polyfit(locs, log_pks, 1);

% Extract D_theta from the slope
%D_theta = -fit_params(1);
% Estimate uncertainty in D_theta (standard error of the slope) :)
%cov_matrix = inv(S.R) * inv(S.R)'; % Covariance matrix of the fit
%D_theta_error = sqrt(cov_matrix(1, 1)); % Uncertainty in the slope

% Generate the fitted curve
%fitted_curve = exp(polyval(fit_params, lag_times)); 

% Display Results
%disp(['Estimated D_theta WITH ERRORS: ', num2str(D_theta), ' ± ', num2str(D_theta_error)]);

%% --- Fit Exponential Decay to Peaks Using Only log_pks_error ---
%-----------------------------------------------------------------

% Compute weights for the linear fit (based on log_pks_error)
weights = 1 ./ (log_pks_error.^2); % Weights are inversely proportional to variance

% Perform weighted linear fit using locs and log_pks
X = [ones(size(locs)), locs]; % Design matrix for linear regression
W = diag(weights); % Weight matrix
beta = (X' * W * X) \ (X' * W * log_pks); % Weighted least squares solution

% Extract D_theta and intercept
D_theta = -beta(2); % Slope is -D_theta
intercept = beta(1); % Intercept of the fit

% Estimate uncertainty in D_theta
cov_matrix = inv(X' * W * X); % Covariance matrix of the fit parameters
D_theta_error = sqrt(cov_matrix(2, 2)); % Uncertainty in the slope (D_theta)


% Show Results of D_theta
disp(['Estimated D_theta (using log_pks_error only): ', num2str(D_theta), ' ± ', num2str(D_theta_error)]);


%% ---Plot Graph---
%--------------------
figure;
plot(lag_times, C, '-o', 'DisplayName', 'C(\Delta t)');
hold on;
plot(lag_times, fitted_curve, '--', 'DisplayName', 'Fit: exp(-D_\theta \Delta t)');
xlabel('\Delta t [s]');
ylabel('C(\Delta t) [unitless]');
%title('Fit to Peaks of Orientation Correlation Function');
legend show;
grid on;


