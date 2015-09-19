function [ smooth_wave ] = db_rectify_waveform(song_waveform, sampling_freq, smoothing_window )
%db_rectify_waveform Input a waveform, sampling frequency, and
%smoothing_window in msecs, return a smoothed and rectified waveform
%   Input your

%rectify waveform
rec_wave = abs(song_waveform);

%create window for smoothing (running average)
window_length = sampling_freq/1000*smoothing_window;
window = ones(1,window_length)./window_length;

smooth_wave = filter(window,1,rec_wave);


end

