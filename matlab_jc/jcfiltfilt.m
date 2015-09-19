function y=jcfiltfilt(x)
[b,a] = butter(4,2/1000,'low');
y=filtfilt(b,a,x);