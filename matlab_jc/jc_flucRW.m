function [k]=jc_flucRW
fs = 32000; % cbin sampling rate
t = 0:1/fs:0.1; %3200pts - 100ms

    a=rand;
    if a>0.5
        signs=1;
    else
        signs=-1;
    end   
k(1)=signs*rand*20;
for i=2:length(t)
    a=rand;
    if a>0.5
        signs=1;
    else
        signs=-1;
    end    
    k(i)=k(i-1)+signs*rand;
end

%x = vco(k,[2940 3160],fs); % 2400 and 3600 are the minimum and maximum frequency values that the pitch oscillates between
