fs=44100


figure

subplot(4,2,1)
taif=xcorr(smaif,smaif,5000);
plot(-5000:5000,taif)
title('original')

subplot(4,2,2) 
tp50m33=xcorr(smaif,smp50m33,5000);
plot(-5000:5000,tp50m33);
title('+50, then - 33');
%axis([0 2 1 10000]);


subplot(4,2,3);
tp2=xcorr(smaif,smp2,5000);
plot(-5000:5000,tp2);

%plot(smp2)%,spectp2,tp2,fp2]=evspect(datp2(:,1),fs);
%title('+ 2');
%axis([0 2 1 10000]);

subplot(4,2,5);
tp10=xcorr(smaif,smp10,5000);
plot(-5000:5000,tp10)
%axis([0 2 1 10000]);

subplot(4,2,7);
tp50=xcorr(smaif,smp50,5000);
plot(-5000:5000,tp50);
%axis([0 2 1 10000]);

subplot(4,2,4);
tm2=xcorr(smaif,smm2,5000);
plot(-5000:5000,tm2);
%axis([0 2 1 10000]);

subplot(4,2,6);
tm10=xcorr(smaif,smm10,5000);
plot(-5000:5000,tm10);
%axis([0 2 1 10000]);

subplot(4,2,8);
tm50=xcorr(smaif,smm50,5000);
plot(-5000:5000,tm50);
%axis([0 2 1 10000]);

