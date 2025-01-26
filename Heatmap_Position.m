% --- Select and Load the CSV File ---
[file, path] = uigetfile({'*.csv', 'CSV Files (*.csv)'}, 'Select the CSV File with Particle Data');
if isequal(file, 0)
    disp('No file selected. Exiting...');
    return;
end
full_path = fullfile(path, file);
data = readtable(full_path);

% --- Extract Particle Data ---
% Ensure your CSV column names match these or adjust accordingly
frame = data.Frame; % Extract frame numbers
x_data = data{:, contains(data.Properties.VariableNames, 'X_Particle')}; % All X_Particle columns from CSV
y_data = data{:, contains(data.Properties.VariableNames, 'Y_Particle')}; % All Y_Particle columns from CSV

% --- Prepare Positions for Heatmap ---
% Reshape data into one-dimensional arrays, ignoring NaN values
x_positions = x_data(:); % Convert X data to a single column
y_positions = y_data(:); % Convert Y data to a single column

% Remove NaN values (if particles are not detected in some frames)
valid_indices = ~isnan(x_positions) & ~isnan(y_positions);
x_positions = x_positions(valid_indices);
y_positions = y_positions(valid_indices);

% --- Create Heatmap ---
x_bins = 0:10:2000; % Adjust grid resolution as needed (10 px bins)
y_bins = 0:10:2000; % Adjust grid resolution as needed (10 px bins)
heatmap_data = histcounts2(y_positions, x_positions, y_bins, x_bins);

% --- Plot Heatmap ---
figure;
imagesc(x_bins, y_bins, heatmap_data); % Plot heatmap
axis equal;
set(gca, 'YDir', 'normal'); % Correct y-axis direction
colormap('hot'); % 'Hot' colormap
colorbar; % Show color scale

% Set the color range to focus on specific values 
%caxis([0 ,50]); % Adjust the limits of the color scale 


%title('Particle Accumulation Heatmap');
xlabel('x [px]');
ylabel('y [px]');

% Set black background
set(gca, 'Color', 'k');
