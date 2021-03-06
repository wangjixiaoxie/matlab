clust=[1 2];
numclust=length(clust);
figure;
plothist=1;
binsize=.005
COLORS='rkmb'
%stimleng=6.18


load /cobain4/twarren4/g100o55/stim/stim.mat
rawsong=corpshiftednormg{1};
stimleng=length(rawsong)/fs;

ax(1)=subplot(numclust*2+1,1,1)
[sm,sp,t,f]=evsmooth(rawsong,44100,0.01);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
axis off;
box off;

if(plothist==0)
    for i=1:numclust
        cur_clst=clust(i);
        ax(i+1)=subplot(numclust*2+1,1,2*i:2*i+1);
        plotrasters3(stimf(cur_clst,1).rast);
        axis([0 stimleng 0 stimf(cur_clst,1).cnt+2])
        if(i==numclust)
            xlabel('Time (s)','Fontsize',16);
            ylabel('Trial','Fontsize', 16,'Color','k')
        end
        if (i~=numclust)
            set(gca,'XTick',[])
            set(gca,'xcolor','w')
            set(gca,'ycolor','k')
        end
    end
    linkaxes(ax(1:numclust+1),'x')
else
    axcount=1;
    for i=1:numclust
        cur_clst=clust(i);
        ax(i+1)=subplot(numclust*2+1,1,2*i:2*i+1);
        edges=0:binsize:stimleng;
        [axt(axcount:axcount+1),h1,h2]=plotrasters3(stimf(cur_clst,1).rast,edges,(1/binsize)*stimf(cur_clst,1).meanhist,stimleng);
        set(h1,'Color','r')
        set(h2,'Color','k')
        for j=1:numclust 
            cur_clst2=clust(j)
            maxval=(1/binsize)*max(stimf(cur_clst2,1).meanhist);
            maxplot=max(maxval)
        end
        axes(axt(axcount))
            axis([0 stimleng 0 stimf(cur_clst,1).cnt+2])
        if(i==numclust)
            XLabel('Time (s)','Fontsize',16);
            YLabel('Trial','Fontsize', 16,'Color','k')
        end
        axes(axt(axcount+1))
        axis([0 stimleng 0 ceil(maxplot)])
        if(i==numclust)
            YLabel('Spike Rate (Hz)','Fontsize', 16,'Color','k')
            
        end
        if (i~=numclust)
            set(gca,'XTick',[])
            set(gca,'xcolor','w')
            set(gca,'ycolor','k')
        end
    axcount=axcount+2
    end
end
 linkaxes(ax(1:numclust+1),'x')   
 linkaxes(axt,'x');


            
         