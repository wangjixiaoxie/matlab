%Make histograms
%outfilefilter={'D:\tim\g36g67stim\3667.wav' 'D:\tim\g36g67stim\3667_rev.wav' 'D:\tim\g36g67stim\3667_m5.wav' 'D:\tim\g36g67stim\3667_p5.wav'...
 %               'D:\tim\g36g67stim\3667l.wav' 'D:\tim\g36g67stim\3667q.wav' 'D:\tim\g36g67stim\3667_p10.wav' 'D:\tim\g36g67stim\3667_p10tl.wav' '3667_m10.wav' '3667_m10tl.wav'}

clear edges            
clustnum=1;                
colorvec={'r' 'y'};            
%each condition
for(ii=[1 2 3 4 5 6 7 8])
    %This is only because of error in stimulus
    
    
    %50ms bins
    clear vec;
    
    %vec=spk1ind(spksin1clust,3);
    indices=find((spkind(spksin{clustnum},3))==ii);
    vec=spkind(spksin{clustnum}(indices),2);
    vec2=length(unique(vec))
    if(ii==1)
    edges{ii}=88203:1600:length(data);
    else
    edges{ii}=1:1600:length(data);
    end
    histout{ii}=histc(spkind(spksin{clustnum}(indices),1),edges{ii});
    histout{ii}=histout{ii}/vec2
end
    
figure
%subplot(2,1,1)
%[sm,sp,t,f]=evsmooth(data(:,song_chan),fs,100);
%imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
%axis([2 10 1 7000])
%subplot(2,1,2)
plot(edges{1}/fs-2.7563,histout{1},'r')
hold on;
%plot(edges/fs,histout{2},'b')
%plot(edges/fs,histout{4},'y')
%plot(edges/fs,histout{5},'m')

plot(edges{2}/fs,histout{7},'b')
%plot(edges{2}/fs,histout{},'g')
plot(edges{2}/fs,histout{9},'k')
box off
xlabel='time'

%plot(edges/fs,histout{6},'y')
%plot(edges/fs,histout{4},'w')
%axis([2 10 0 6])

legend('bos','p10','m10')
legend boxoff