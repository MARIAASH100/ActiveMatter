%% initial parametars
%filename = '..\FRISBEE_18';
%vid_filename= strcat(filename,'.mp4'); % *replace* with your video name
N_particles=14; % *replace* with number of particles

% --- Select Video File using Dialog ---
[file, path] = uigetfile({'*.mp4;*.avi;*.mov', 'Video Files (*.mp4, *.avi, *.mov)'}, 'Select a Video File');
if isequal(file, 0)
    disp('No video selected. Exiting...');
    return;
end

% --- Construct Full Path to the Selected Video ---
full_path = fullfile(path, file);
%v = VideoReader(full_path); % Read video



%% load video from file
v = VideoReader(full_path); % Read video
%v = VideoReader(vid_filename); % read video
NF=v.NumFrames; % the number of frames
t_start=v.CurrentTime; % initial starting timestamp

X_mean_t=zeros(1,NF); % define arrays for X and Y positions of the center of mass
Y_mean_t=zeros(1,NF);
N_detected_t=zeros(1,NF); % define array for the number of particles identified in each frame
X_all=zeros(N_particles,NF); % define arrays for X and Y positions of all particles
Y_all=zeros(N_particles,NF);

%% display a single frame to find particle radius
frame = readFrame(v); % this command reads the video frames in order and advances the internal clock 'CurrentTime'
figure(); hold on; box on;
imshow(frame) % this command displays the current frame as an image
v.CurrentTime=t_start; % we rewind the video to its beginning

%% choose parameters
R_min=10; % R_min, R_max in pixels for Hough transform (circle detection limits)
R_max=30;
Pixel_to_cm=29; % conversion for using cm units; *replace* with your value
fps=30; % conversion of time to seconds; *replace* with your fps value

%% loop over frames and identify circles
v.CurrentTime=t_start;

for ii=1:NF % we loop over all frames in the video
    cur_frame = readFrame(v); % read current frameth
    [~,masked_frame]=frame_isolate_green_channel(cur_frame); % we apply a custom mask to filter only the green stickers
    masked_frame_green=masked_frame(:,:,2); % after applying the mask we take only the green channel for the RGB representation
    masked_frame_green_rescaled=masked_frame_green.*(255/double(max(max(masked_frame_green)))); % we rescale the green values between 0 and 255 for maximal precision
    [centers, radii] = imfindcircles(masked_frame_green_rescaled,[R_min R_max],'ObjectPolarity','bright',Sensitivity=0.92); % this function identifies circles using a circular Hough transform
    N_detected=length(radii);
    if N_detected>1
        remove=[];
        for jj=1:N_detected-1 % if two identified particles are closer than the sum of their radii - we mark one of them for removal (misidentification)
            for kk=jj+1:N_detected
                r_ij=((centers(jj,1)-centers(kk,1))^2+(centers(jj,2)-centers(kk,2))^2)^0.5;
                if r_ij<(radii(jj)+radii(kk))
                    remove=[remove,kk];
                end
            end
        end
        remove=unique(remove); % we consider only unique values for removal
        radii(remove)=[]; % we remove the flagged particles
        centers(remove,:)=[];
        N_detected=N_detected-length(remove); % this is the final number of identified articles
    end
    if N_detected>0
        N_detected_t(ii)=N_detected; % fill the relevant spot in the arrays we defined
        X_mean_t(ii)=mean(centers(:,1),1);
        Y_mean_t(ii)=mean(centers(:,2),1);
        X_all(1:N_detected,ii)=centers(:,1);
        Y_all(1:N_detected,ii)=centers(:,2);
    end

    
%         % this next part displays every frame along with the identified
%         % particles. you can use it to make sure your image processing is working
%         % well

%         figure(); hold on; box on;
%         imshow(cur_frame)
%         viscircles(centers, radii,'EdgeColor','b');
%         title(strcat('frame=',string(ii),';  time=',string(v.CurrentTime),';  particles_detected=',string(length(centers))))

end

%% exporting to csv:
% Prepare data for export
frame_indices = (1:NF)'; % Frame numbers
time_stamps = (frame_indices - 1) / fps; % Time stamps in seconds

% Initialize particle coordinate columns
particle_data = cell(N_particles * 2, 1);
for p = 1:N_particles
    particle_data{2*p-1} = X_all(p, :)'; % X coordinate for particle p
    particle_data{2*p} = Y_all(p, :)'; % Y coordinate for particle p
end

% Generate flat variable names
particle_var_names = {};
for p = 1:N_particles
    particle_var_names = [particle_var_names, {['X_Particle_' num2str(p)], ['Y_Particle_' num2str(p)]}];
end

% Combine all variable names into a flat array
all_var_names = [{'Frame', 'Time_s'}, particle_var_names, {'ParticlesDetected', 'X_CenterOfMass', 'Y_CenterOfMass'}];

% Combine all data into a table
export_data = table(frame_indices, time_stamps, ...
    particle_data{:}, N_detected_t', X_mean_t', Y_mean_t', ...
    'VariableNames', all_var_names);

% Write to CSV file
csv_filename = fullfile(path, 'particle_tracking_data.csv'); % Save in the same directory as the video
writetable(export_data, csv_filename);

disp(['Data exported to ', csv_filename]);


