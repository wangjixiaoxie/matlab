function [ db_waveform ] = db_waveform_to_decibels(waveform, sampling_freq, smoothing_window )
%db_wavform_to_decibels Converts amplitude in the given waveform to
%decibels by finding
%   Detailed explanation goes here

%calculate window length (convert from msecs)
window = sampling_freq/1000*smoothing_window;

%sliding rms calculation
for i = 1:length(waveform)-window
    rms_wave(i) = rms(waveform(i:i+window));
end

min_rms = min(rms_wave);

db_waveform = 20*log10(rms_wave./min_rms);


end

