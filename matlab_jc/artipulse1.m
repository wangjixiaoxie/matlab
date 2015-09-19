function [artipuls,articheck]=artipulse1(harmonic)
fs=32000;
for i=1:30
    
    %Random pitch shifting so that residuals has to work
    b=rand;
    if b<0.5
        c=-1;
    else
        c=1;
    end
    d=0.5+0.5*b;
    c=c*d;
    
    
    phase=rand*0.12;
    frequency=160;
    t=0:1/fs:0.12;
    articheck=c*0.3*sin(frequency*pi*(t+phase))+0.3*c*sin(40*pi*t-phase);
    for n=1:1000;
        articheck(n)=1*rand;
    end
    for n=length(articheck)-1000:length(articheck)
        articheck(n)=1*rand;
    end
    %0.33*sin(20*pi*(t+0.05))+0.33*sin(60*pi*t)+0.33*sin(1000*pi*t);
    artipuls(i).pitches=vco(articheck,[harmonic*1200 harmonic*3500],fs);

    %[ifdgram,sonogram]=ifdv(x,f3162s,1024,1020,1,1,3,5,5);

end