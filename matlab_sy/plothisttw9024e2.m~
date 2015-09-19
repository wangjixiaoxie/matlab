figure
clear maxtrial;
slength=8.331;
load '/cobain3/twarren3/w90o24e2/songs.mat'


imagevecar{1}=imagevec;
imagevecar{2}=altimagevec;
snums=1:8
%subplotind=[3:7 9:10 12]


%Following lines are a bit of a hack to make sure image of shifted stimulus
%turns out all right

zeroind=find(imagevec==0);
oneind=find(imagevec==1);
negind=find(imagevec==-1);


imagevecar{1}=imagevec;
imagevecar{1}(negind)=1;
imagevecar{1}(zeroind)=1;
imagevecar{1}(oneind)=0;


zeroind=find(altimagevec==0);
oneind=find(altimagevec==1);
negind=find(altimagevec==-1);

imagevecar{2}=altimagevec;
imagevecar{2}(negind)=1;
imagevecar{2}(zeroind)=1;
imagevecar{2}(oneind)=0;


imagevals=[8 11];
subplotvals=[3:7 9:10 12]
structnames={'bos' 'rev' 'm10' 'p10' 'jig' 'jigc' 'alt' 'altc'}
rawwav='/cobain3/twarren/w90o24e2/9024.wav';
splen=length(subplotvals)+4;



ax(1:2)=subplot(splen,1,1:2);

    

    [sm,sp,t,f]=evsmooth(corpshiftednormg{1},44100,0.01);
    imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
    box off;
    hold on;

    %for ii=1:length(xlist)
     %   plot([xlist(ii) xlist(ii)],[0 10000],'w');hold on;
    %end
    %set(gca,'Xtick',xlist);
    %set(gca,'Xticklabel',1:length(xlist));
    box off; axis off;

for kk = 1:length(snums)%length(stimf)
	
    
    
    
    
    ax(subplotvals(kk))=subplot(splen,1,subplotvals(kk));	
        xv=[1:size(stimf(snums(kk)).hist,2)]*histbinwid;
		yv=stimf(snums(kk)).hist./stimf(snums(kk)).cnt./histbinwid;
		maxtrial(snums(kk))=max(yv);
        
        
        plot(xv,yv,[clr(1),'-']);
        
        box off;
        %axis off;
        ylabel(structnames{kk},'Fontsize',16);
     
          if (kk~=length(snums))
            
            set(gca,'xcolor','w')
          
            set(gca,'YTickLabel',[]);
            set(gca,'Ytick',[]); 
            set(gca,'XTickLabel',[]);
            set(gca,'Xtick',[]);
            
          end
          if(kk==length(snums))
              
            xlabel('Time(s)', 'Fontsize', 16);
            set(gca,'xcolor','k')
          end
end



%This is to scale the data correctly.
for kk=1:length(snums)
    ax(subplotvals(kk))=subplot(splen,1,subplotvals(kk));
    hold on;
    maxval=max(maxtrial);
    plot(0:.1:slength, maxval/2, 'r--')
    
    axis([0 slength 0 maxval]);
end

for ii=1:length(imagevals)
        
        imaget=[0:.001:length(imagevecar{1})/1000-.001]
        ax(imagevals(ii))=subplot(splen,1,imagevals(ii));
        %clim=[[0 0 0] ;[1 1 1]]
        imagesc(imaget,f,imagevecar{ii})
        colormap('gray')
        %ax(imagevals(ii))=subplot(splen,1,imagevals(ii));
        %xlabel('Time (s)','Fontsize',16);
        %set(gca,'YTickLabel','Fontsize',14);
        %set(gca,'Xtick',[]);
        box off; axis off;
        %set(gca,'Ytick',[]);
        %ylabel('jigstim','Fontsize',14)
    end

    %box off;
    
    linkaxes(ax,'x');