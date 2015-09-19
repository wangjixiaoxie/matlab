figure
colormap(hot)
path1='g21g81_0180307_0711.13964.cbin';
%path2='g21g81_010207_0808.345.cbin';
cd /doyale1/twarren/g21g81/stim25
[dat1,fs]=evsoundin('',path1,'obs0');
rec1=readrecf(path1)
%[dat2,fs]=evsoundin('',path2,'obs0');
%rec2=readrecf(path2)
%plot(0:1/fs:(length(data(:,2))-1)/fs,data(:,2))
figure;ax(1)=subplot(8,1,1:3)
[sm,sp,t,f]=evsmooth(dat2,fs,100);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
colormap(hot)
box off; axis on;ylabel('Frequency (Hz)')
ax(2)=subplot(8,1,4)
x=[rec2.ttimes/1000;rec2.ttimes/1000]
y=[10;20]
plot(x,y,'r','Linewidth',3)
box off; axis off;
title(path1)
ax(3)=subplot(8,1,5:7)
[sm,sp,t,f]=evsmooth(dat1,fs,100);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
box off;
title(path2)
ax(4)=subplot(8,1,8)
x=[rec1.ttimes(1)/1000;rec1.ttimes(1)/1000];
plot(x,y,'r','Linewidth',3)
box off; axis off;
linkaxes(ax,'x')
axis([0 6 10 20])