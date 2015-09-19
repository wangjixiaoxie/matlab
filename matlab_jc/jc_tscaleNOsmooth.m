function [a,xax]=jc_tscaleNOsmooth(PC,LATENT)
a=zeros(1,25);
for i=1:25
    xax(i)=(1/(i/2))*length(PC)/8;
end
for i=1:5
    
    hh=jc_fft(xcorr(PC(:,i)));
    for k=1:25
        a(k)=a(k)+LATENT(i)*hh(k);
    end
end