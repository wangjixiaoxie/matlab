function [ filter_song ] = db_bandpass_waveform(waveform, sampling_freq, freq_low, freq_high )
%db_bandpass_waveform Uses a bandpass filter on the raw waveform.
%   Input waveform, sampling frequency, and low/high cutoff frequencies.
%   Uses a butterworth filter

%checks to see if freq_high is too high, if so changes it
if freq_high >= .5*sampling_freq-500
    disp('current freq_high is too high')
    freq_high = .5*sampling_freq-1000;
    disp(['reset to ' freq_high ' Hz'])
end

%creates 8 pole butterworth filter
[b, a] = butter(8,[freq_low*2/sampling_freq, freq_high*2/sampling_freq]);
filter_song = filtfilt(b, a, waveform);




end

