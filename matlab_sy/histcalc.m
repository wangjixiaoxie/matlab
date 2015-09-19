%Make histograms
%outfilefilter={'D:\tim\g36g67stim\3667.wav' 'D:\tim\g36g67stim\3667_rev.wav' 'D:\tim\g36g67stim\3667_m5.wav' 'D:\tim\g36g67stim\3667_p5.wav'...
 %               'D:\tim\g36g67stim\3667l.wav' 'D:\tim\g36g67stim\3667q.wav' 'D:\tim\g36g67stim\3667_p10.wav' 'D:\tim\g36g67stim\3667_p10tl.wav' '3667_m10.wav' '3667_m10tl.wav'}

            
clustnum=1;                
colorvec={'r' 'y'};            
%each condition
for(ii=[1:10])
    %50ms bins
    clear vec;
    
    %vec=spk1ind(spksin1clust,3);
    indices=find((spkind(spksin{clustnum},3))==ii);
    vec=spkind(spksin{clustnum}(indices),2);
    
    vec2=length(unique(vec))
    edges=1:1600:length(data);
    histout{ii}=histc(spkind(spksin{clustnum}(indices),1),1:1600:length(data));
    histout{ii}=histout{ii}/vec2
end
    
figure
ax(1)=subplot(2,1,1)
[sm,sp,t,f]=evsmooth(data(:,song_chan),fs,100);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
%axis([2 10 1 7000])
ax(2)=subplot(2,1,2)
plot(edges/fs,histout{1},'r')
hold on;
plot(edges/fs,histout{2},'g')
plot(edges/fs,histout{3},'b')
%plot(edges/fs,histout{6},'g')

%plot(edges/fs,histout{8},'g')
%plot(edges/fs,histout{3},'m')
box off
xlabel='time'

%plot(edges/fs,histout{6},'y')
%plot(edges/fs,histout{4},'w')
%axis([2 10 0 6])

legend('bos','m10tl','p10tl','Color','k')%,'p10','m10')
legend boxoff

linkaxes(ax(1:2),'x');
gca(ax(1));
axis([0 7 100 10000])