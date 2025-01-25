% Load the video
%v = VideoReader('D:\Programs in D\Matlab active matter\Part1VideosSingle\c5.mp4');
% --- Select Video File using Dialog ---
[file, path] = uigetfile({'*.mp4;*.avi;*.mov', 'Video Files (*.mp4, *.avi, *.mov)'}, 'Select a Video File');
if isequal(file, 0)
    disp('No video selected. Exiting...');
    return; 
end

% --- Construct Full Path to the Selected Video ---
full_path = fullfile(path, file);
v = VideoReader(full_path); % Read video

% Read the first frame
frame = read(v, 1);

% Display the frame
figure;
imshow(frame);
title('Select the top-left and bottom-right corners of the region of interest');

% Use ginput to select two points
[x, y] = ginput(2);

% Convert points to integers (row and column indices)
row_start = round(y(1));
row_end = round(y(2));
col_start = round(x(1));
col_end = round(x(2));

% Show the selected region
roi = frame(row_start:row_end, col_start:col_end, :); % Region Of Interest 
figure;
imshow(roi);
title('Selected Region of Interest');

% Save the ROI coordinates for further use
fprintf('Row range: %d to %d\n', row_start, row_end);
fprintf('Column range: %d to %d\n', col_start, col_end);
