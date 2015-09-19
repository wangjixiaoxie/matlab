function [ amplitude_envelope, time, fs, hilbert_envelope, normal_song] = db_song_amplitude_envelope( filename, window_size, raw )
%db_song_amplitude_envelope Input a cbin file and it will calculate the
%amplitude envelope (either max or average) with a window size (default
%8 ms).
%   Detailed explanation goes here

if nargin < 2
    window_size = 8;
    raw = 0;
end

%Gets the wavefore, time, and sampling rate from the song
[normal_song, time, fs, raw_song, bandpass_raw] = db_normalize_song_waveform(filename);

%Calculates the number of bins to get a window of time window_size
window = (window_size*fs)/1000;

%uses a hilbert transformation to calculate amplitude envelope
if raw == 0
    hilbert_song = hilbert(normal_song);
else
    hilbert_song = hilbert(bandpass_raw);
end

hilbert_envelope = sqrt(hilbert_song .* conj(hilbert_song));

%uses runningaverage to smooth the hilbert_envelope
amplitude_envelope = smooth(hilbert_envelope,window);
amplitude_envelope = amplitude_envelope';

% %%% calculates the rate of change of the hilbert_envelope
% diff_amp_env = zeros(length(amplitude_envelope),1);
% for i = 2:length(amplitude_envelope)-1
%     diff_amp_env(i) = (amplitude_envelope(i+1) - amplitude_envelope(i-1))./(2/fs);
% end
% diff_amp_env = diff_amp_env./max(diff_amp_env);
% 
% hilbert_diff = hilbert(diff_amp_env);
% hilbert_diff_env = sqrt(hilbert_diff .* conj(hilbert_diff));
% 
% der_amp_env = runningaverage(hilbert_diff_env,window);
% der_amp_env = der_amp_env';
%         
    
    


end

