%Make histograms
%outfilefilter={'D:\tim\g36g67stim\3667.wav' 'D:\tim\g36g67stim\3667_rev.wav' 'D:\tim\g36g67stim\3667_m5.wav' 'D:\tim\g36g67stim\3667_p5.wav'...
 %               'D:\tim\g36g67stim\3667l.wav' 'D:\tim\g36g67stim\3667q.wav' 'D:\tim\g36g67stim\3667_p10.wav' 'D:\tim\g36g67stim\3667_p10tl.wav' '3667_m10.wav' '3667_m10tl.wav'}


bintimes=[2 7]; 
clustnum=1;                
colorvec={'r' 'y'};            
numplots=10
edges=0:100:2000
%each condition

%make a simple matrix, which calculates for every trial-- how many spikes
%occurred in the bin
%loop through each trial, find what the stimulus is
%add the spike number to a vector of numbers
for stim=1:10
    spikedist{stim}=[];
end


for (ifn=1:max(spkind(:,2)))
    ind=find(spkind(:,2)==ifn);
    stim=spkind(ind(1),3);
    spksinbin=calcspikesinbin(spkind(ind,1),bintimes,fs)
    spikedist{stim}=[spikedist{stim};spksinbin];
end





for(stim=[1:numplots])
    %50ms bins
    clear vec;
    
    %vec=spk1ind(spksin1clust,3);
    
    
    
    vec2=length(spikedist{stim})
    
    
    
    histspikedist{stim}=histc(spikedist{stim},edges);
    histspikedist{stim}=histspikedist{stim}/vec2;
end
    
figure
ax(1)=subplot(2,1,1)
[sm,sp,t,f]=evsmooth(data(:,song_chan),fs,100);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
%axis([2 10 1 7000])
ax(2)=subplot(2,1,2)
a=colormap;
for stim=1:3
    plot(edges,histspikedist{stim},'Color',a(floor(stim/3*64),:),'Linewidth',3)
    hold on;
end
box off
xlabel='NUM SPIKES IN BIN'

%plot(edges/fs,histout{6},'y')
%plot(edges/fs,histout{4},'w')
%axis([2 10 0 6])

%legend(structnames)
legend('bos', 'rev', 'loud')
legend boxoff

%linkaxes(ax(1:2),'x');
%gca(ax(1));
%axis([0 7 100 10000])