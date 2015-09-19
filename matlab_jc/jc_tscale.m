function [a,b]=jc_tscale(PC,LATENT)
j=0;
for i=1:5
    hh=jc_fft(xcorr(PC(:,i)));
    for k=1:length(hh)
        t(k)=i;
    end
    [xf]=lsqcurvefit(@jc_2gauss,[100 3 1 50 5 1 10 9 1],t,hh');
    a(j+1)=(1/(xf(2)*0.5))*25;
    b(j+1)=LATENT(i)*xf(1);
    j=j+1;
    a(j+1)=(1/(xf(5)*0.5))*25;
    b(j+1)=LATENT(i)*xf(4);
    j=j+1;
    a(j+1)=(1/(xf(8)*0.5))*25;
    b(j+1)=LATENT(i)*xf(7);
    j=j+1;
end