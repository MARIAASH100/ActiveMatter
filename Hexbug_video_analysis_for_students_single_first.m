%cd 'C:\Users\USER\Documents\Active_Matter\students\Anna Maria\single_particle'
%filename = '..\c5';
%vid_filename= strcat(filename,'.mp4'); % *replace* with your video name
%v = VideoReader(vid_filename); % read video

% --- Select Video File using Dialog ---
[file, path] = uigetfile({'*.mp4;*.avi;*.mov', 'Video Files (*.mp4, *.avi, *.mov)'}, 'Select a Video File');
if isequal(file, 0)
    disp('No video selected. Exiting...');
    return;
end

% --- Construct Full Path to the Selected Video ---
full_path = fullfile(path, file);
v = VideoReader(full_path); % Read video

NF=v.NumFrames; % the number of frames
t_start=v.CurrentTime; % initial starting timestamp



x_t=zeros(1,NF);  % define arrays for X and Y positions of the particle
y_t=zeros(1,NF);
t = zeros(1,NF);

%%
frame = readFrame(v); % this command reads the video frames in order and advances the internal clock 'CurrentTime'
figure(); hold on; box on;
imshow(frame) % this command displays the current frame as an image
v.CurrentTime=t_start; % we rewind the video to its beginning

%%
min_r=30; % R_min, R_max in pixels for Hough transform (circle detection limits)
max_r=200;
%%

v.CurrentTime=t_start; % we rewind the video to its beginning

for ii=1:NF % we loop over all frames in the video
    t(ii) = v.CurrentTime;

    frame = readFrame(v);
    frame_gray=rgb2gray(frame);
    frame_gray_norm=(frame_gray-min(min(frame_gray)))*(255/((max(max(frame_gray))-min(min(frame_gray)))));

    [centers, radii] = imfindcircles(frame_gray_norm,[min_r max_r],'ObjectPolarity','dark'); %Hough transform assuming dark circles
%      figure();
%      imshow(frame_gray_norm)
% 
%      viscircles(centers, radii,'EdgeColor','b');
    try
        x_t(ii)=centers(1);
        y_t(ii)=centers(2);
        disp(ii)
    catch
        continue
    end
end

%%
figure;
plot(x_t, y_t);
grid on;


%% Export Data to CSV
output_data = [(1:NF)', x_t', y_t']; % Combine frame numbers, X, and Y into a matrix
csv_filename = 'trajectory_data_with_frame_numbers.csv'; % Specify the output file name
writematrix(output_data, csv_filename); % Write the matrix to a CSV file

disp(['Trajectory data saved to ', csv_filename]);