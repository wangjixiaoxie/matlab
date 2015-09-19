function plotLMANxcorr(x1,x2)
% e.g. x1=ACSF(3).residHILB2(1900:2500,21:50);
% e.g. x2=INA(3).residHILB2(2000:2600,1:30);

% 1. produce the autocorrelation of the "LMAN-dependent variability"

% Power spectral density method
xcx1=mean(jc_xcorr(x1)');
xcx2=mean(jc_xcorr(x2)');
xcz1=fft(xcx1);
xcz2=fft(xcx2);
xcmag1=abs(xcz1);
xcmag2=abs(xcz2);
xcpha1=angle(xcz1);
xcpha2=angle(xcz2);
xcmag=(xcmag1-xcmag2);
xcp=xcpha2;
xcre=xcmag.*cos(xcp);
xcim=xcmag.*sin(xcp);
xcdat=complex(xcre,xcim);
xcgg=real(ifft(xcdat));
hold on;plot(xcgg./max(xcgg),'k')