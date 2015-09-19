inarray=[];
figure
clustnum=1;
numstim=6
binsize=.005
stim2plot= [1  6 5 4 2 3 ];
labels={'bos' 'sn' 'p1.2' 'm1.2'  'm2'  'rev'};
drawx=0;
fs=32000;
%rawsong=data(34.8*32000:40.96*32000,1);
%make a vector which can be imaged,
%which is 1-d and goes from start points to stop points in time.

%respvec=conscresp{clustnum}./(1/binsize);

stimleng=stmpclen;

%make a vector of meanhists

for i=1:length(stim2plot)
    stm=stim2plot(i);
    inarray=[inarray;stimf(clustnum,stm).meanhist];
end
%this is to plot vector with pitchshift pattern
%plothistfxn(inarray/binsize,binsize,labels,rawsong,thresh/binsize,stimleng,respvec)

    plothistfxn(inarray/binsize,binsize,labels,rawsong,44053,stimleng)
    if(drawx) 
        subplot(numstim+1,1,1)
        hold on;
        ind=find(randvecfin==1);
        xvals=xlist(ind);
        xvals2=xlist(ind+1);
        yarray=2500;
        drawxrasters(xvals',xvals2',yarray,'r')
    end    
        