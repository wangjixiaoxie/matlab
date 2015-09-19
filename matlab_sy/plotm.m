NCLUST = 1%length(stimf);

ii=1

    clr=['brkgmcy'];
    hh=figure;hold on;
    for kk = 1:1%length(stimf)
        figure;
        %xv=[1:size(stimf(kk).rast,2)]*binwid;
        %yv=[1:stimf(kk).cnt];
        %imagesc(xv,yv,stimf(kk).rast);syn;
         box off;
        
       
        
        
        edges=0:.05:stimleng;
        %edges2=[0 startsong*32000 endsong*32000]
        for (stim=1:1)
        
            trials=stimf(stim).cnt;
            for trialnum=1:trials
                ind2=find(stimf(stim).spkarray(:,2)==trialnum);
                if(isempty(ind2))
                    stimf(stim).histdist(trialnum,length(edges))=0;
                else
                    stimf(stim).histdist(trialnum,:)=(histc(stimf(stim).spkarray(ind2,1),edges))
        %histdist3{stim}(trialnum,:)=(histc(spkind(spksin{clustnum}(ind(ind2)),1),edges2))
                end
            end
            end
    end   
 
        ax1 = gca;
        
        stimf(stim).meanhist=mean(stimf(stim).histdist)
        
         
        [axt,h1,h2]=plotrasters3(stimf(kk).rast,edges,20*stimf(stim).meanhist,stimleng);
        axes(axt(1))
        axis([0 stimleng 0 trials+1]) 
        ylabel('Trial','Fontsize',16,'Color','r');
        
        axes(axt(2))
        axis([0 stimleng 0 100])
        ylabel('Firing Rate (Hz)','Fontsize',16);
        xlabel('Time (s)','Fontsize',16);
        set(gca,'YTick',[0 50 100])
        
        
        
   
        
        
        
        
        
        
        
   