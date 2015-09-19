function [x]=jc_fluca
fs = 32000; % cbin sampling rate
t = 0:1/fs:0.1; %3200pts - 100ms
for i=1:length(t)
    k(i)=0;
end
x = vco(k,[4900 5100],fs); % 2400 and 3600 are the minimum and maximum frequency values that the pitch oscillates between

%y=noise*rand(size(x));  % add in some noise at each point in time

%figure; spectrogram(x,128,124,128,32000, 'yaxis')

%Gaussian windows (sonogram)
%sigma=2;
%OVERLAP=1010;
%N=1024;
%t=-N/2+1:N/2;
%sigma=(sigma/1000)*fs;
%w=exp(-(t/sigma).^2);
%q=specgram(x,1024,[],w,OVERLAP)+eps;
%sonogram=abs(flipdim(q,1));