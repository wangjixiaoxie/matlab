low_tone=rand(32000*0.05,1)*2-1;
wavwrite(white,32000,16,'wn_16khz50ms.wav');

%%
Fs = 32000;      %# Samples per second
toneFreq = 1000;  %# Tone frequency, in Hertz
nSeconds = 0.1;   %# Duration of the sound
y = sin(linspace(0, nSeconds*toneFreq*2*pi, round(nSeconds*Fs)));
% y2 = 0.6*sin(linspace(0, nSeconds*2*toneFreq*2*pi, round(nSeconds*Fs))); % 2nd harmonic. 0.6x intensity
% y=y1+y2;
% sound(y, Fs);  %# Play sound at sampling rate Fs


%%
wavwrite(y, Fs, 16, 'tone_1000Hz_100ms.wav');  %# Save as an 8-bit, 1 kHz signal