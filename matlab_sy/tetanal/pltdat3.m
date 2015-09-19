


tet_chans=[3:6];song_chan=1;

flstoplot=[15]
titlelist={ 'trial 15'}
chanstoplot=[5 6]



for ii=1:length(flstoplot)
 fn=fnm(flstoplot(ii)) 
  [data,fs,spkindnew,spkampnew]=tetanal(fn{1},-2000,song_chan,tet_chans);
   datar{ii}=data;
end




chanstoplot=[5 6]
a=[1:10]

figure

ax=subplot(3,1,1)
[sm,sp,t,f]=evsmooth(datar{1}(:,song_chan),fs,100);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
axis([2 10 1 7000])



for ii=1:length(flstoplot)
    
    
    indices=find(spk1ind(spksin2clust,2)==flstoplot(ii));
    spks2toplot=spk1ind(spksin2clust(indices),1);
    ampstoplot2c3=spk1amp(spksin2clust(indices),3);
    ampstoplot2c4=spk1amp(spksin2clust(indices),4);
    
    
    
    indices=find(spk1ind(spksin1clust,2)==flstoplot(ii));
    spks1toplot=spk1ind(spksin1clust(indices),1);
    
    
    
    ampstoplot1c3=spk1amp(spksin1clust(indices),3);
    ampstoplot1c4=spk1amp(spksin1clust(indices),4);
    
    
    numtoplot1=length(spks1toplot);
    numtoplot2=length(spks2toplot);
    
    data=datar{ii};
    ax(2*ii)=subplot(3,1,2*ii)
    plot((1:length(data(:,5)))/fs,data(:,5),'k','LineWidth',1)
    hold on;
    x1=[(spks1toplot'/fs); (spks1toplot'/fs)]
    y1=[-ampstoplot1c3'.*ones(1,numtoplot1)-1500;-ampstoplot1c3'.*ones(1,numtoplot1)-500];
    x2=[spks2toplot'/fs; spks2toplot'/fs]
    y2=[-ampstoplot2c3'.*ones(1,numtoplot2)-1500;-ampstoplot2c3'.*ones(1,numtoplot2)-500];
    
    
    plot(x1,y1,'r','LineWidth',2)
    plot(x2,y2,'b','LineWidth',2)
    
    %title({titlelist{ii} 'chan3'} )
    
  axis([6.85 6.9 -10000 10000])
   ax(2*ii+1)=subplot(3,1,2*ii+1)
   plot((1:length(data(:,6)))/fs,data(:,6),'k','LineWidth',1)
   hold on;
   x1=[spks1toplot'/fs; spks1toplot'/fs]
    y1=[-ampstoplot1c4'.*ones(1,numtoplot1)-1500;-ampstoplot1c4'.*ones(1,numtoplot1)-500];
    x2=[spks2toplot'/fs; spks2toplot'/fs]
    y2=[-ampstoplot2c4'.*ones(1,numtoplot2)-1500;-ampstoplot2c4'.*ones(1,numtoplot2)-500];
   
    plot(x1,y1,'r','LineWidth',2)
    plot(x2,y2,'b','LineWidth',2)
   
   
   axis([6.85 6.9 -10000 10000])
end





    