function [xax]=summary_plot2(FFs,batch,avgpitch)
fs=32000;
[vals,trigs]=triglabel(batch,'a',1,1,1,0);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);


[avna,t,f]=get_avn(batch,'a',0,0.2,'','','obs0'); 

for i=1:length(avgpitch)
    xax(i)=((i+128)/8)/1000;
end
figure;hold on;
imagesc(t,f,log(avna));syn;ylim([0,1e4]);
plot(xax,avgpitch,'g')
plot(toff/1000,mean(FFs),'*')
plot(0.06,FFs,'*')
xlim([0 0.18]);
ylim([1 4000]);