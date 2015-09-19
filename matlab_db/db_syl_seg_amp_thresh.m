function [onset, offset ] = db_syl_seg_amp_thresh( waveform, fs, min_duration, min_gap, threshold, varargin )
%db_syl_seg_amp_thresh Summary of this function goes here
%   Detailed explanation goes here

%checks the number of arguments, if option smoothing_window is not
%specified, it will be set to 2 msecs
if nargin < 6
    smoothing_window = 2;
else
    smoothing_window = varargin{1};
end

%converts min_duration & min_gap to fs length
min_duration = fs/1000*min_duration;
min_gap = fs/1000*min_gap;

%bandpass filters waveform
filtered_song = db_bandpass_waveform(waveform, fs, 500, 10000);

%normalize song
normalized_song = filtered_song./max(abs(filtered_song));

%converts filtered_song in a.u. to decibels (min decibel is set as lowest
%value of song with duration smoothing_window).
%final_wv = db_waveform_to_decibels(normalized_song, fs, smoothing_window);

%rms waveform instead of decibel
smoothing_window = fs/1000*smoothing_window;
final_wv = zeros(length(normalized_song)-smoothing_window,1);
for i = 1:length(normalized_song)-smoothing_window
    final_wv(i) = rms(normalized_song(i:i+smoothing_window));
end

%find onsets
j=1;
for i = min_gap+2:length(final_wv)-min_duration
    if all(final_wv(i-min_gap:i-1) < threshold) &&...
            all(final_wv(i:i+min_duration) > threshold)
        onset(j) = i;
        j = j+1;
    end
end

%find offsets
k=1;
for i = min_duration+2:length(final_wv)-min_gap
    if all(final_wv(i-min_duration:i-1) > threshold) && ...
            all(final_wv(i:i+min_gap) < threshold)
        offset(k) = i;
        k = k+1;
    end
end



end

