function [xax]=adhocplotter926(batch,avgpitch)


[avna,t,f]=get_avn(batch,'a',0,0.2,'','','obs0'); 

for i=1:length(avgpitch)
    xax(i)=((i+128)/8)/1000;
end
figure;hold on;
imagesc(t,f,log(avna));syn;ylim([0,1e4]);
plot(xax,avgpitch,'g')
xlim([0 0.18]);
ylim([1 4000]);