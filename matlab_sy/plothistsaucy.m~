figure

xind=0:.02:(length(stimf(1).meanhist)-1)/50
subplot(5,1,1)
ax2(5)=plot(0:1/32000:(length(dat(:,1))-1)/32000, dat(:,1))


ax2(1)=subplot(5,1,2);
plot(xind,stimf(1).meanhist)
box off;
title('bos');
ax2(2)=subplot(5,1,3);
plot(xind,stimf(2).meanhist)
box off;
title('rev');
ax2(3)=subplot(5,1,4);
plot(xind,stimf(3).meanhist)
box off;
title('m10');
ax2(4)=subplot(5,1,5);
plot(xind,stimf(4).meanhist)
box off;
title('p10')
linkaxes(ax2,'x')
linkaxes(ax2(1:4))
