%ONLY FIX FOR G7980.
%histcalc7 was created 1/5/05, deals with different data file sizes




%Make histograms
%outfilefilter={'D:\tim\g36g67stim\3667.wav' 'D:\tim\g36g67stim\3667_rev.wav' 'D:\tim\g36g67stim\3667_m5.wav' 'D:\tim\g36g67stim\3667_p5.wav'...
 %               'D:\tim\g36g67stim\3667l.wav' 'D:\tim\g36g67stim\3667q.wav' 'D:\tim\g36g67stim\3667_p10.wav' 'D:\tim\g36g67stim\3667_p10tl.wav' '3667_m10.wav' '3667_m10tl.wav'}

clear histdist;
clear histdist3;
clear norm_meanallspikes;
clear error_allspikes

bintimes=[2 7]; 
clustnum=1;                
colorvec={'r' 'b' 'k' 'g'};            
numplots=5;
ratio=0;
%startsong=5.7;
%endsong=6.1;
%datalength=8;
%each condition

%make a simple matrix, which calculates for every trial-- how many spikes
%occurred in the bin
%loop through each trial, find what the stimulus is
%add the spike number to a vector of numbers
for stim=1:5
    histdist{stim}=[];
end


for (stim=1:numplots)
    
    ind=find(spkind(spksin{clustnum},3)==stim);
    if(max(spkind(spksin{clustnum}(ind),1))>400000)
        datalength(stim)=12;
        startsong(stim)=6;
        endsong(stim)=12;
    else
        datalength(stim)=10;
        startsong(stim)=4;
        endsong(stim)=10;
    end
    edges{stim}=0:1600:datalength(stim)*32000;
    edges2{stim}=[0 startsong(stim)*32000 endsong(stim)*32000]

    trials=unique(spkind(spksin{clustnum}(ind),2));
    for trialnum=1:length(trials)
        ind2=find(spkind(spksin{clustnum}(ind),2)==trials(trialnum));
        histdist{stim}(trialnum,:)=(histc(spkind(spksin{clustnum}(ind(ind2)),1),edges{stim}))
        histdist3{stim}(trialnum,:)=(histc(spkind(spksin{clustnum}(ind(ind2)),1),edges2{stim}))
    end
   
   end

%now find mean and variance of this across trials






for(stim=[1:numplots])
    if(isempty(histdist{stim})==0)
    
    
        numtrials=length(histdist{stim}(:,1))
        norm_meanallspikes{stim,:}=mean(histdist{stim})
        error_allspikes{stim,:}=std(histdist{stim})/sqrt(numtrials)
        rat{stim}=(histdist3{stim}(:,2)/(endsong(stim)-startsong(stim)))./(histdist{stim}(:,1)/2);
        if(ratio==0)
            rat{stim}=(histdist3{stim}(:,2)/(endsong(stim)-startsong(stim)));
        end
        
        meanspikerat(stim)=mean(rat{stim});
        errspikerat(stim)=std(rat{stim})/sqrt(numtrials);

    end
end    
figure
ax(1)=subplot(2,1,1)
[sm,sp,t,f]=evsmooth(data(:,song_chan),fs,100);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
axis([2 10 1 10000])
ax(2)=subplot(2,1,2)
a=colormap;

figure
stim=[1 2 3 4 5 ]
for i=1:4
    
    if(length(edges{stim(i)})>201)
         plot(edges{stim(i)}(41:241)/fs-2,norm_meanallspikes{stim(i),1}(41:241),'Color',colorvec{i},'Linewidth',2)
    else
        
        plot(edges{stim(i)}/fs,norm_meanallspikes{stim(i),:},'Color',colorvec{i},'Linewidth',1)
    end
    hold on;
    %plot(edges/fs,norm_meanallspikes(stim,:)+error_allspikes(stim,:),'--','Color',colorvec{stim},'Linewidth',1)
    %plot(edges/fs,norm_meanallspikes(stim,:)-error_allspikes(stim,:),'--','Color',colorvec{stim},'Linewidth',1)
end
axis([0 12 0 4])
box off
xlabel='NUM SPIKES IN BIN'

%plot(edges/fs,histout{6},'y')
%plot(edges/fs,histout{4},'w')
%axis([2 10 0 6])

%legend(structnames)
legend('bos', 'rev', 'm10', 'p10','lou')
legend boxoff
figure

errorbar(meanspikerat([1 2 3 4 5]) ,errspikerat([1 2 3 4 5 ]),'+','Linewidth',3)
box off 
set(gca,'xtick',[1:9]);
%set(gca,'xticklabel',['bos';'rev';'l  ';'m10';'m10';'m5 ';'p5 ';'m1t';'p1t';'qui']);
set(gca,'xticklabel',['bos';'rev';'m10';'p10';'lou'],'Fontsize',14);
set(gca,'Ylabel','Ratio of song SR to presong SR') 


%linkaxes(ax(1:2),'x');
%gca(ax(1));
%axis([0 7 100 10000])