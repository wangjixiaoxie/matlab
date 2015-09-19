function [artipuls,articheck]=artipulse2(harmonic)
fs=32000;
for k=1:30

blips=20;

%Random sign (+/-)
b=rand;
if b<0.5
   c=-1;
else
   c=1;
end
    
width=1400;
%Random disalignment of the slow (LMAN) component with the center of the note
slow_shift=round(c*rand*width); %The number of points it gets shifted to the right of the center of the "note"
medium_shift=round(c*rand*width);
for i=1:blips
    d=rand;
    if d<0.5
        e=-1;
    else
        e=1;
    end
    fast_shift(i)=round(e*rand*width);
end


%Width of the slow (LMAN) component
sd=0.6; %Std of the slow Gaussian component
coeff_slow=1/(2*sd*sd); %Gaussian coefficient

sdm=0.05;
coeff_medium=1/(2*sdm*sdm);

%Width of the fast component
sdf=0.005;
coeff_fast=1/(2*sdf*sdf);

%Weights
slow_weight=0.4;
fast_weight=0.03;
medium_weight=0.1+rand*0.1;

t=-1:10/fs:1;
for i=(width+1):length(t)-(width+1);
    Slow(i)=slow_weight*c*exp(-coeff_slow*t(i-slow_shift)*t(i-slow_shift));
    Medium(i)=medium_weight*c*exp(-coeff_medium*t(i-medium_shift)*t(i-medium_shift));
    %0.1*exp(-400*t(i-1200)*t(i-1200))+0.1*exp(-400*t(i+400)*t(i+400))+0.1*
    %exp(-200*t(i+100)*t(i+100))+0.1*c*exp(-20*t(i)*t(i));
end

%Add in a blip
for i=(width+1):length(t)-(width+1);
    for j=1:blips
        fast(j,i)=fast_weight*exp(-coeff_fast*t(i-fast_shift(j))*t(i-fast_shift(j)));
    end
end

total_fast=fast(1,:);
for i=2:blips
    total_fast=total_fast+fast(i,:);
end


articheck=total_fast+Medium+Slow;
articheck=articheck(3200-width:3200+width);
%figure; plot(articheck); ylim([mean(articheck)-7*std(articheck) mean(articheck)+7*std(articheck)]);
artipuls(k).pitches=vco(articheck,[harmonic*1200 harmonic*3500],fs);
end