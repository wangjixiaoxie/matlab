function lagt=lineup(datnorm,datshift,fs);
%lagt=lineup(datnorm,datshift,fs);
[sm1]=evsmooth(datnorm,fs,0);
[sm2]=evsmooth(datshift,fs,0);

[c,l]=xcorr(sm1,sm2,5000);                       
[y,i]=max(c);
lagt=l(i)/fs;
return;
