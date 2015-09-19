fs=44100


figure



subplot(4,2,1)
[smaif,spectaif,taif,faif]=evspect(dataif,fs,0);
title('original');
axis([0 2 1 10000]);

subplot(4,2,2)
[smp50m33,spectp50m33,tp50m33,fp50m33]=evspect(datp50m33(:,1),fs,-.0176);
title('+50, then - 33');
axis([0 2 1 10000]);


subplot(4,2,3);
[smp2,spectp2,tp2,fp2]=evspect(datp2(:,1),fs,-.0166);
title('+ 2');
axis([0 2 1 10000]);

subplot(4,2,5);
[smp10,spectp10,tp10,fp10]=evspect(datp10(:,1),fs,-.0089);
title('+ 10');
axis([0 2 1 10000]);

subplot(4,2,7);
[smp50,spectp50,tp50,fp50]=evspect(datp50(:,1),fs,.0088);
title('+ 50');
axis([0 2 1 10000]);

subplot(4,2,4);
[smm2,spectm2,tm2,fm2]=evspect(datm2(:,1),fs,-.0207);
title('- 2');
axis([0 2 1 10000]);

subplot(4,2,6);
[smm10,spectm10,tm10,fm10]=evspect(datm10(:,1),fs,-.0236);
title('- 10');
axis([0 2 1 10000]);

subplot(4,2,8);
[smm50,spectm50,tm50,fm50]=evspect(datm50(:,1),fs,-.0314);
title('- 50');
axis([0 2 1 10000]);

