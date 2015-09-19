% Requires some input vector d5 that contains the chunk of the syllable of
% which you want to take the fft.  Code is taken directly from the matlab
% help menu describing the function fft.
samp=44100; % assumes sampling rate is 44100
window=441; % size of the area you want the fft over (e.g. 441=10ms)
window_shift=220; %shift between each window (e.g. 220=5ms)
window_left=1-window_shift;
for i=1:8
    window_left=window_left+window_shift;
    window_right=window_left+window;
    chunk=d5(window_left:window_right); % chunk of data you're taking the fft of
    NFFT = 2^nextpow2(window); % Next power of 2 from length of y
    Y = fft(chunk,NFFT)/window;
    f = samp/2*linspace(0,1,NFFT/2); %fs=44100

    % Plot single-sided amplitude spectrum.
    subplot(2,4,i); plot(f,2*abs(Y(1:NFFT/2))) 
end