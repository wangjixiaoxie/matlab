function [ final_onset, final_offset, run_time ] = db_syllable_segmentation( waveform, fs, compress, volume_thresh, min_duration, max_duration, L )
%db_syllable_segmentation Summary of this function goes here
%   Detailed explanation goes here
tic
if nargin < 3
%     b = 0.2;
%     k = .6;
    compress = 0;
    volume_thresh = 0.01;
    min_duration = 10;
    max_duration = 500;
    L = .01;
elseif nargin < 4
    volume_thresh = 0.01;
    min_duration = 10;
    max_duration = 500;
    L = .01;
elseif nargin < 5
    min_duration = 10;
    max_duration = 500;
    L = .01;
elseif nargin < 6
    max_duration = 500;
    L = .01;
elseif nargin < 7
    L = .01;
end

%% Bandpass filter song (50,10000 Hz)
fprintf('\nFiltering and normalizing song...')

filter_song = db_bandpass_waveform(waveform, fs, 500, 10000);
time = linspace(0,length(filter_song)/fs,length(filter_song))';
%Normalize sob = 0.01;ng so range is -1 to 1
% normalizing_factor = max(abs(filter_song));
% normalized_song = filter_song./normalizing_factor;
fprintf('done!\n')
    
%% Calculate rms (64 points ~ 2 ms for 32kHz) with 50% overlap
fprintf('\nCalculating rectified waveform...')
window = 32;
overlap = .75;

% rms_waveform = zeros(length(1:round(overlap*window):length(normalized_song)-window),1);
% time_rms = time(1:round(overlap*window):length(normalized_song)-window);
%j = 1;

[rms_waveform, indices] = db_rms(filter_song,window,round(window*overlap),1);

time_rms = time(indices);
rms_waveform = rms_waveform';

%normalize rms_waveform
rms_waveform = rms_waveform./max(rms_waveform);

%dynamic compression of waveform
if compress == 1;
    for i = 1:length(rms_waveform)
        if rms_waveform(i) > 0.5
            rms_waveform(i) = rms_waveform(i)*0.5;
        end
    end
    
    %normalize waveform
    rms_waveform = rms_waveform./max(rms_waveform);
end

j = 1;
min_window =  fs/1000*100;
for i = 1:min_window:length(rms_waveform)-min_window
    min_wave(j) = mean(rms_waveform(i:i+min_window));
    j = j+1;
end

rms_waveform = rms_waveform-min(min_wave);

% for i = 1:round(overlap*window):length(normalized_song)-window
%     rms_waveform(j) = db_rms(normalized_song(i:i+window));
%     j = j+1;
% end

fprintf('done!\n')
%% Calculate forward and backward vectors (rate of change of envelope)
%calculates a,b for two line segments radiating from each point (best fit
%line). Number of points needed to calculate a,b is determined by adding
%the euclidean distance between two adjacent points until a threshold is
%reached (L). This allows for few points to be necessary when there is a sudden
%change in amplitude, but many points when there is little change.
fprintf('\nCalculating change in amplitude envelope (may take a while)...')
back_vec = zeros(2,length(rms_waveform));
forw_vec = back_vec;
angle = zeros(length(rms_waveform),1);
threshold_cross  = angle;

for i = 2:length(rms_waveform)-1
    
    %calculates number of points BEFORE current point to fit line
    Epsilon_back = 0;
    j_back = i;
    counter_back = 0;
    while Epsilon_back < L && j_back > 1
        Epsilon_back = Epsilon_back + norm([rms_waveform(j_back+1);time_rms(j_back+1)]-[rms_waveform(j_back);time_rms(j_back)]);
        j_back = j_back-1;
        counter_back = counter_back + 1;
    end
    
    %calculates number of points AFTER current point to fit line
    Epsilon_for = 0;
    j_for = i;
    counter_forward = 0;
    while Epsilon_for < L && length(rms_waveform)-j_for >= 1
        Epsilon_for = Epsilon_for + norm([rms_waveform(j_for+1);time_rms(j_for+1)]-[rms_waveform(j_for);time_rms(j_for)]);
        j_for = j_for+1;
        counter_forward = counter_forward + 1;
    end
    
    %least squares fit for back vector
    back = db_least_squares_lin_fit(time_rms(j_back:i)-time_rms(i), rms_waveform(j_back:i)-rms_waveform(i));
    %least square vector centered at origin
    back_vec(:,i) = [-1,back(1)*-1+back(2)];
    
    %least squares fit for forward vector
    forw = db_least_squares_lin_fit(time_rms(i:j_for)-time_rms(i), rms_waveform(i:j_for)-rms_waveform(i));
    %least square vector centered at origin
    forw_vec(:,i) = [1,forw(1)*1+forw(2)];
    
    %display(i)
end

fprintf('done!\n')
%% calculates the angle between forward vector and backward vector
fprintf('\nCalculating threshold crossing...')

for i = 2:length(rms_waveform)
    angle(i) = 2*pi-db_angle_2d_vectors(forw_vec(:,i),back_vec(:,i));
    %display(i)
end

for i = 2:length(rms_waveform)
    threshold_cross(i) = abs(angle(i)-pi);
end

%since this angle is highly sensitive to even very small changes in
%amplitude, we set all threshold_crosses to zero if it is not above some
%very minimal amplitude threshold.
for i = 1:length(threshold_cross)-window
    if rms_waveform(i) < volume_thresh + min(rms_waveform)
        threshold_cross(i) = 0;
    end
end


fprintf('done!\n')
%% onsets calculation
fprintf('\nCalculating onsets...')
%threshold for change in angle between forward and backward vectors
angle_thresh = .5; %max is pi
%minimum duration for a syllable
min_duration_fs = round(fs/1000*min_duration*length(time_rms)/length(time));

%finds all threshold_crosses above the angle thresh
above_angle_thresh = find(threshold_cross > angle_thresh);
%finds the time differences between elements above
time_between_thresh_cross = diff(above_angle_thresh);
%sets up onsets vectors (will be in indices of rms_waveform)
onsets = zeros(length(find(time_between_thresh_cross > min_duration_fs))+1,1);
%first onset is first time threshold_cross crosses angle_thresh
onsets(1) = above_angle_thresh(1);
%next onsets are after min_duration gap in time_between_thresh_cross
onsets(2:end) = above_angle_thresh(find(time_between_thresh_cross > min_duration_fs)+1);

onset_times = time_rms(onsets);

fprintf('Done!\n')
%% offsets calculation
fprintf('\nCalculating offsets...')
%min threshold for change in angle between forward and backward vectors
%min_angle_thresh = 1.5;
%minimum gap duration
%min_gap_fs = fs/1000*min_gap;

%calculates offsets, by looking between onset(i) and onset(i+1) when
%threshold_cross goes below min_angle_thresh for min_gap time.
offsets = zeros(length(onsets),1);
if length(onsets) > 1
    for i = 1:length(onsets)-1
        above_min_angle = find(threshold_cross(onsets(i):onsets(i+1)) > angle_thresh);
        offsets(i) = onsets(i) + above_min_angle(end-1);
    end
end

%for last syllable in onsets
above_min_angle = find(threshold_cross(onsets(end):end) > angle_thresh);
offsets(end) = onsets(end) + above_min_angle(end-1);

offset_times = time_rms(offsets);
fprintf('done!\n')
%% gets rid of any 'syllable' that is too short or too long (min_duration, max_duration)
fprintf('\nFiltering out noise...')
j = 1;
for i = 1:length(onsets)
    if (offset_times(i)-onset_times(i))*1000 > min_duration && (offset_times(i)-onset_times(i))*1000 < max_duration
        offset_filter1(j) = offset_times(i);
        onset_filter1(j) = onset_times(i);
        j = j+1;
    end
end

offset_filter1 = offset_filter1';
onset_filter1 = onset_filter1';

%% gets rid of any 'syllable' where the area underneath the rms_waveform is very small (< 1)
j = 1;
for i = 1:length(onset_filter1)
    if trapz(rms_waveform(find(time_rms == onset_filter1(i)):find(time_rms == offset_filter1(i)))) > 1
        
        final_onset(j) = onset_filter1(i);
        final_offset(j) = offset_filter1(i);
        j = j+1;
    end
end

final_onset = final_onset';
final_offset = final_offset';

fprintf('Done!\n')
run_time = toc;
end

