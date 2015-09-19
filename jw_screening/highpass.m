function y=highpass(x,cutoff,sr)
%usage: y=highpass(x,cutoff,sr);
%high-pass  cutoff filter
%5th order butterworth non-causal; b=nth order 

[b,a]=butter(5, cutoff/sr,'high');
y=filtfilt(b,a,x);



%function y=lowpass(x,Fc,Fs)
%%usage: y=lowpass(x);
%%high-pass 200 Hz cutoff filter
%%5th order butterworth non-causal; b=nth order 
%%hard-wired for 16000 hz sampled data
%
%[b,a]=butter(5, Fc/Fs);
%y=filtfilt(b,a,x);
