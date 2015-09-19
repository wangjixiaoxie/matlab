%Make histograms
%outfilefilter={'D:\tim\g36g67stim\3667.wav' 'D:\tim\g36g67stim\3667_rev.wav' 'D:\tim\g36g67stim\3667_m5.wav' 'D:\tim\g36g67stim\3667_p5.wav'...
 %               'D:\tim\g36g67stim\3667l.wav' 'D:\tim\g36g67stim\3667q.wav' 'D:\tim\g36g67stim\3667_p10.wav' 'D:\tim\g36g67stim\3667_p10tl.wav' '3667_m10.wav' '3667_m10tl.wav'}

clear histdist;
clear norm_meanallspikes;
clear error_allspikes

bintimes=[2 7]; 
clustnum=2;                
colorvec={'r' 'b' 'k' 'g'};            
numplots=4;
startsong=4;
endsong=8;
datalength=9;
%each condition

%make a simple matrix, which calculates for every trial-- how many spikes
%occurred in the bin
%loop through each trial, find what the stimulus is
%add the spike number to a vector of numbers
for stim=1:10
    histdist{stim}=[];
end

edges=0:1600:datalength*32000;
edges2=[0 startsong*32000 endsong*32000]
for (stim=1:numplots)
    ind=find(spkind(spksin{clustnum},3)==stim);
    trials=unique(spkind(spksin{clustnum}(ind),2));
    for trialnum=1:length(trials)
        ind2=find(spkind(spksin{clustnum}(ind),2)==trials(trialnum));
        histdist{stim}(trialnum,:)=(histc(spkind(spksin{clustnum}(ind(ind2)),1),edges))
        histdist2{stim}(trialnum,:)=(histc(spkind(spksin{clustnum}(ind(ind2)),1),edges2))
    end
        
end

%now find mean and variance of this across trials






for(stim=[1:numplots])
    
    numtrials=length(histdist{stim}(:,1))
    norm_meanallspikes(stim,:)=mean(histdist{stim})
    error_allspikes(stim,:)=std(histdist{stim})/sqrt(numtrials)
    rat{stim}=(histdist2{stim}(:,2)/5)%./(histdist2{stim}(:,1)/2);
    meanspikerat(stim)=mean(rat{stim});
    errspikerat(stim)=std(rat{stim})/sqrt(numtrials);
end
    
figure
ax(1)=subplot(2,1,1)
[sm,sp,t,f]=evsmooth(data(:,song_chan),fs,100);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
axis([2 9 1 10000])
ax(2)=subplot(2,1,2)
a=colormap;
for stim=[1 2 3 4]
    plot(edges/fs,norm_meanallspikes(stim,:),'Color',colorvec{stim},'Linewidth',3)
    hold on;
    %plot(edges/fs,norm_meanallspikes(stim,:)+error_allspikes(stim,:),'--','Color',colorvec{stim},'Linewidth',1)
    %plot(edges/fs,norm_meanallspikes(stim,:)-error_allspikes(stim,:),'--','Color',colorvec{stim},'Linewidth',1)
end
axis([2 9 0 2])
box off
q=Xlabel('Time (s)','Fontsize', 15);
q=Ylabel('Average Spikes in 50ms bin', 'Fontsize',15);
%plot(edges/fs,histout{6},'y')
%plot(edges/fs,histout{4},'w')
%axis([2 10 0 6])

%legend(structnames)
legend('bos', 'rev', 'p10', 'm10')
legend boxoff
figure

errorbar(meanspikerat([1:4]) ,errspikerat([1:4]),'+','Linewidth',3)
box off 
set(gca,'xtick',[1:4]);
%set(gca,'xticklabel',['bos';'rev';'l  ';'m10';'p10';'m5 ';'p5 ';'m1t';'p1t';'qui']);
set(gca,'xticklabel',['bos';'rev';'p10';'m10'], 'FontSize',15);
p=Ylabel('Average Spike Rate during song (Hz)','FontSize',15) 


%linkaxes(ax(1:2),'x');
%gca(ax(1));
%axis([0 7 100 10000])