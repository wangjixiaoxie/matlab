function [TH]=get_thresh2(data,song_chan,tet_chans,chn2plot,TH,fs,fn);
%data is scaled to actual units for this thresh
%plot the data,when the user hits return, prompt a ginput statment

ptstoplt=1:64000;


figure;

yvec=data(ptstoplt,chn2plot);
xvec=0:1/fs:(length(yvec)-1)/fs;
xvec2=0:1/100:(length(yvec)-1)/fs;
subplot(2,1,1)
plot(data(:,1));
title(fn);
subplot(2,1,2)
plot(xvec,yvec);
hold on;
%[ans]=input('press return when zoomed in');
if(TH(1)==0)

    [x,y]=ginput();

    mn_dat=mean(data(floor(x(1)*fs):ceil(x(2)*fs),chn2plot));
    st_dat=std(data(floor(x(1)*fs):ceil(x(2)*fs),chn2plot));


    TH(1)=mn_dat+6*st_dat;
    TH(2)=mn_dat+12*st_dat;

    
end
TH
thvec1=TH(1)*ones(length(xvec2),1);
thvec2=TH(2)*ones(length(xvec2),1);
plot(xvec2,thvec1','r','Linewidth',2);
plot(xvec2,thvec2','r','Linewidth',2);
ans=input('Threshold okay?')
if(isempty(ans))
   
else
    [x,y]=ginput();
    TH=y;
    
end  