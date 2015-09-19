[dat,fs]=readevtaf(fn,'0r');
[sm,sp,t,f]=evsmooth(dat,fs,0.01);      
imagesc(t,f,log(abs(sp)));set(gca,'YD','n');
m=colormap('gray');colormap(m(end:-1:1,:));hold on;
[trig,fs]=readevtaf(fn,'1r');plot([1:length(trig)]/fs,trig*3000,'r');
hold off;
