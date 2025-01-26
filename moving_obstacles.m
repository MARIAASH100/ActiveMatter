
% --- Select Video File using Dialog ---
[file, path] = uigetfile({'*.mp4;*.avi;*.mov', 'Video Files (*.mp4, *.avi, *.mov)'}, 'Select a Video File');
if isequal(file, 0)
    disp('No video selected. Exiting...');
    return;
end

% --- Construct Full Path to the Selected Video ---
full_path = fullfile(path, file);
video = VideoReader(full_path); % Read video

% Initialize variables
frameJump = 300; % Number of frames to skip
cm_per_pixel = 1 / 29; % Conversion factor from pixels to centimeters
frameRate = video.FrameRate; % Frames per second

% Read the first frame as the reference frame
video.CurrentTime = 0; % Start from the first frame
frame = readFrame(video);

% Show the first frame and then select a reference point (with your hands)
figure;
imshow(frame);
title('Select the initial reference point');
[x1, y1] = ginput(1); % You select the initial reference point
x1 = round(x1);
y1 = round(y1);

% Initialize lists for distances and time
distances = [];
timeStamps = [];

% Loop through frames with the defined jump
frameNumber = 1; % Start at the first frame
while hasFrame(video)
    if mod(frameNumber, frameJump) == 1
        % Read and display the current frame
        frame = readFrame(video);
        imshow(frame);
        title(sprintf('Frame %d - Select a point', frameNumber));
        
        % Select a point on the frame(rest)
        [xi, yi] = ginput(1);
        xi = round(xi);
        yi = round(yi);
        
        % Calculate the distance from the reference point
        distancePixels = sqrt((xi - x1)^2 + (yi - y1)^2);
        distanceCm = distancePixels * cm_per_pixel;
        
        % Store the distance and corresponding time
        distances = [distances, distanceCm];
        timeStamps = [timeStamps, frameNumber / frameRate];
    else
        % Skip the frame
        readFrame(video);
    end
    frameNumber = frameNumber + 1;
end

% Plot d(t) in cm
figure;
plot(timeStamps, distances, '-o');
xlabel('t [s]');
ylabel('Distance [cm]');
%title('Distance over Time');
grid on;



% Plot d(t) in cm with log scale for the y-axis
%figure;
%semilogy(timeStamps, distances, '-o'); % Log scale for the y-axis
%xlabel('t [s]');
%ylabel('Distance [cm]');
%title('Distance over Time (Log Scale)');
%grid on;

disp('Processing complete!');