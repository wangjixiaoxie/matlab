function [vlscomb,sumdata]=plotwnlearningfigv4(vlscomb)
%first arm up vls file
%modified to emphasize pitchshift
% figstoplot=[ 1  3 4];


% figstoplot=[1 3 4 10];
figstoplot=[4  ]
threshdays=[0 1 1 2 2 8]
threshvals=[2.44 2.44 2.47 2.47 2.52 2.52]
% if(~exist('vlscomb'))
% 
% %oct 29 and oct30
% pth{1}='/oriole2/bk41bk14/518comb'
% cvl{1}='batch18.catch'
% 
% pth{2}='/oriole2/bk41bk14/astim519'
% cvl{2}='batch20.catch'
% 
% pth{3}='/oriole2/bk41bk14/astim519'
% cvl{3}='batch22.catch'
% 
% 
% 
% %aug22
% pth{4}='/oriole2/bk41bk14/a-bstim523'
% cvl{4}='batch24.catch'
% 
% pth{5}='/oriole2/bk41bk14/wnon516'
% cvl{5}='batch17.catch'
% 
% pth{6}='/oriole2/bk41bk14/astim519'
% cvl{6}='batch19.catch'
% 
% pth{7}='/oriole2/bk41bk14/a-bstim523'
% cvl{7}='batch26.catch'
% 
% pth{8}='/oriole2/bk41bk14/astim519'
% cvl{8}='batch21.catch'
% 
% pth{9}='/oriole2/bk41bk14/a-bstim523'
% cvl{9}='batch25.catch'
% 
% %pitchinds
% inds{1}=[1:6 8:9]
% 
% %syntax
% pth{13}='/oriole2/bk41bk14/temptest2'
% cvl{13}='batch03.catch.keep'
% 
% pth{14}='/oriole2/bk41bk14/wnonsyn1'
% cvl{14}='batch06.catch'
% 
% pth{15}='/oriole2/bk41bk14/wnonsyn1'
% cvl{15}='batch08.catch'
% 
% 
% pth{16}='/oriole2/bk41bk14/temptest3'
% cvl{16}='batch04.catch.keep'
% 
% pth{10}='/oriole2/bk41bk14/wnonsyn1'
% cvl{10}='batch05.catch'
% 
% pth{11}='/oriole2/bk41bk14/wnonsyn1'
% cvl{11}='batch07.catch'
% 
% pth{12}='/oriole2/bk41bk14/wnonsyn1'
% cvl{12}='batch10.catch'
% %syninds
% inds{2}=10:16
% 
% % pth{4}='/oriole/bk48w74/BASELINE4/temptest/files01'
% % cvl{4}='batch.rand.rand.keep'
% % 
% % pth{5}='/oriole/bk48w74/BASELINE4/temptest/files02'
% % cvl{5}='batchcomb'
% % pth{6}='/oriole/bk48w74/BASELINE4/temptest/files03'
% % cvl{6}='batch.rand.keep'
% % 
% % pth{7}='/oriole/bk48w74/BASELINE4/temptest/files04'
% % cvl{7}='batch.rand.keep'
% % 
% % pth{8}='/oriole/bk48w74/BASELINE4/temptest/files05'
% % cvl{8}='batch.rand.keep'
% % pth{9}='/oriole/bk48w74/BASELINE4/temptest/files06'
% % cvl{9}='batch.rand.keep'
% % 
% % pth{10}='/oriole/bk48w74/BASELINE4/temptest/files07'
% % cvl{10}='batch.rand.keep'
% % 
% 
% %notes day 1 thresh varied from 2.42 to 2.46, plot as 2.44
% 
% 
% % pth{9}='/oriole2/bk57w35evren/bk57w35/amponaug19'
% % cvl{9}='batch28.rand'
% % vlscomb=[];
% SOUNDSTR='wav'
% STIMSTR='rig'
% STIMBNDS=[200 -200]
% for analtype=1:2
%     vlscomb{analtype}{1}=[];
%     vlscomb{analtype}{2}=[];
% for ii=1:length(inds{analtype})    
%   crind=inds{analtype}(ii);
%     cmd=['cd ' pth{crind}];
%   eval(cmd);
%   bt=cvl{crind};
%   tbinshft=.074
%   NFFT=512;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
%   fbins=[3200 4000];
%   save BINS_B NFFT fbins tbinshft
%   % frequency analysis just for 'b'
%   load BINS_B
%   NT{1}='a'
%   NT{2}='b'
%   PRENT='';
%   PSTNT='';
%   for ntind=1:2
%     fv{ntind}=findwnote4(bt,NT{ntind},PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
%     
% %      fvst=findwnote9(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0',0,STIMBNDS, SOUNDSTR,STIMSTR);
% %     fvd=findwnote4(bt2,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
%     
%         vls=getvals(fv{ntind},1,'TRIG');
% % for ii=1:length(fv)
% %        if(~isempty(fvst(ii).SOUNDTIME))
% %            vls(ii,4)=fvst(ii).SOUNDTIME;
% %        else
% %            vls(ii,4)=0
% %        end
% %     end
%     
%         vlscomb{analtype}{ntind}=[vlscomb{analtype}{ntind};vls]
%     end
%         sumdata=[];
%     
% %     [vals,trigs]=triglabel(bt,'a',1,1,0,0);
% 
% end
% end
% end
if(ismember(1,figstoplot))
   figure
   sumdata=[];
   subplot(4,4,9:10)
   vlsin=vlscomb{1}{1}
   [unq_floor,mnout]=calcmeanout(vlsin);
   startday=1;
   init_tm=floor(mnout(startday,1))+.5;
   plot(vlsin(:,1)-init_tm,vlsin(:,2)/1000,'o','MarkerSize',1,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor','none');
   hold on;
%    indht=find(vlscomb(:,4)>0&vlscomb(:,4)<105&(vlscomb(:,1)-init_tm>0))
   
%    plot(vlscomb(indht,1)-init_tm,vlscomb(indht,2)/1000,'o','MarkerSize',1,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor','none');
   plot(mnout(:,1)-init_tm,mnout(:,2)/1000,'ko','MarkerSize',4)
    plot(mnout(:,1)-init_tm,mnout(:,2)/1000,'k','Linewidth',1)
    daspect([10 1 1])
    plot(threshdays,threshvals,'r--')
    axis([-1 9 3 3.8])
    daspect([6 1 1])
end

%plot three histograms - init histogram over last baseline day.
%day 3 all
%day 3 hits
if(ismember(2,figstoplot))
%     figure
    subplot(3,4,10:11);
    edges=[6800:43.9394:8250]
    floorvls=floor(vlscomb(:,1));
    indbas=find(floorvls<=unq_floor(2));
    hstbas=histc(vlscomb(indbas,2),edges);
    mnbas=mean(vlscomb(indbas,2))/3000;
    stdbas=std(vlscomb(indbas,2))/3000;
    sumdata.mnbas=mnbas;
    sumdata.stdbas=stdbas;
    sumdata.nbas=length(indbas);
    hstbas=hstbas./length(indbas);
    ind_three=find(floorvls==unq_floor(6));
    hst_three=histc(vlscomb(ind_three,2),edges);
    indht_three=intersect(ind_three,indht);
    hstht=histc(vlscomb(indht_three,2),edges);
    meanht=mean(vlscomb(ind_three,2))/3000;
    stdhit=std(vlscomb(ind_three,2))/3000;
    sumdata.mnhit=meanht;
    sumdata.stdhit=stdhit;
    hstlearn=hst_three./length(ind_three);
    hstht=hstht./length(ind_three);
    edges=edges(1:end);
    hstbas=hstbas(1:end);
    hstlearn=hstlearn(1:end);
    hstht=hstht(1:end);
    
    stairs(edges/3000,hstbas,'k')
    hold on;
    stairs(edges/3000,hstlearn,'b')
%     stairs(edges/3000,hstht,'r')
    plot(mnbas,0.3,'Marker','v','Color','k');
    plot(meanht,0.3,'Marker','v','Color','b');
    unqthresh=unique(threshvals);
    plot([unqthresh;unqthresh],[0 0 0; 0.4 0.4 0.4],'r--')
end  
  
if(ismember(3,figstoplot))
%     figure
    subplot(4,4,13:14);
    edges=[3000:40:3750]
    floorvls=floor(vlsin(:,1));
    indbas=find(floorvls<=unq_floor(1));
    hstbas=histc(vlsin(indbas,2),edges);
    mnbas=mean(vlsin(indbas,2))/1000;
    hstbas=hstbas./length(indbas);
    
    ind_three=find(floorvls==unq_floor(3));
    hst_three=histc(vlsin(ind_three,2),edges);
%     indht_three=intersect(ind_three,indht);
%     hstht=histc(vlsin(indht_three,2),edges);
    meanht=mean(vlsin(ind_three,2))/1000;
    hstlearn=hst_three./length(ind_three);
%     hstht=hstht./length(ind_three);
    edges=edges(1:end-1);
    hstbas=hstbas(1:end-1);
    hstlearn=hstlearn(1:end-1);
%     hstht=hstht(1:end-1);
    
    stairs(edges/1000,hstbas,'k')
    hold on;
    stairs(edges/1000,hstlearn,'b')
%     stairs(edges/3000,hstht,'r')
    plot(mnbas,0.15,'Marker','v','Color','k');
    plot(meanht,0.15,'Marker','v','Color','b');
    axis([3.05 3.7 0 0.3])
    daspect([1 1 1])
end  
    
if(ismember(4,figstoplot))
%     ps.plot_trigs=1;
%     exsong(1).path='/oriole/bk48w74/BASELINE4/temptest/files30'
% %     exsong(1).cvl=cvl{2}
%     exsong(1).fn='bk48w74_301009_170051.397.cbin';
%     exsong(1).bnds=[6 7.9]
%     exsong(1).ax=subplot(3,4,3:4)
%     ps.detect=0;
%     plotcbin(exsong,ps);
    
    exsong(1).path='/oriole2/bk41bk14/screen'
%     exsong(1).cvl=cvl{2}
    exsong(1).fn='bk41bk14_300409_1527.82.cbin';
    exsong(1).bnds=[4.6 7.9]
    exsong(1).ax=subplot(4,4,1:4)
    ps.detect=0;
   [datout]= plotcbin(exsong,ps);
      mxvl=max(abs(datout));
    datoutnrm=datout/(mxvl*1.0001);
    wavwrite(datoutnrm,32000,16,'bk41bk14all.wav')
    ax=gca();
    plotbox([.6 .68 0 10000],'c',ax);
    
    ptsong(1).path='/oriole2/bk41bk14/screen'
    ptsong(1).fn='bk41bk14_300409_1527.82.cbin';
    ptsong(1).bnds=[5.2 5.28]
    ptsong(1).ax=subplot(4,4,5:6);
    [datout]=plotcbin(ptsong,ps);
     mxvl=max(abs(datout));
     datoutnrm=datout/(mxvl*1.0001);
    ax=gca();
    plotbox([.6 .68 0 10000],'c',ax);
     wavwrite(datoutnrm,32000,16,'bk41bk14short.wav')
    
    
    
    synsong(1).path='/oriole2/bk41bk14/wnonsyn1'
    synsong(1).fn='bk41bk14_070509_1104.13516.cbin'
    synsong(1).bnds=[2.9 5.5]
    synsong(1).ax=subplot(4,4,7:8);
    ps.detect=0;
    plotcbin(synsong,ps);
%     exsong(1).path='/oriole/bk48w74/BASELINE4/temptest/files30'
% %     exsong(1).cvl=cvl{2}
%     exsong(1).fn='bk48w74_301009_170051.397.cbin';
%     exsong(1).bnds=[7.4 7.55]
%     exsong(1).ax=subplot(3,4,1)
%     ps.detect=0;
%     plotcbin(exsong,ps);
%     hold on;
%     plot([.086 .086],[0 10000],'r--');
%     hold on;
%      plot([.1020 .1020],[0 10000],'r--');
    
%     exsong(1).stim_path='/oriole/bk48w74/stimexample'
%     exsong(1).stim_mat='stimexample.mat'
end
   

if(ismember(10,figstoplot))
%    figure
   sumdata=[];
   subplot(4,4,15:16)
   [unq_floor,synout]=calcsynout(vlscomb{2});
%    startday=1;
   init_tm=unq_floor(1);
   pcta=synout(:,1)./((synout(:,1)+synout(:,2)));
   pctb=synout(:,2)./((synout(:,1)+synout(:,2)));
   plot(unq_floor-init_tm,pcta*100,'o','MarkerSize',5,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor','none');
   hold on;
   plot(unq_floor-init_tm,pcta*100,'Color',[0.4 0.4 0.4]);
%    hold on;
%    indht=find(vlscomb(:,4)>0&vlscomb(:,4)<105&(vlscomb(:,1)-init_tm>0))
   hold on;
   plot(unq_floor-init_tm,pctb*100,'o','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','none');
   plot(unq_floor-init_tm,pctb*100,'Color','r');
   axis([-1 9 0 100])
   daspect([1 30 1])
%    plot(vlscomb(indht,1)-init_tm,vlscomb(indht,2)/1000,'o','MarkerSize',1,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor','none');
%    plot(mnout(:,1)-init_tm,mnout(:,2)/1000,'ko','MarkerSize',4)
%     plot(mnout(:,1)-init_tm,mnout(:,2)/1000,'k','Linewidth',1)
%     daspect([10 1 1])
%     plot(threshdays,threshvals,'r--')
end



function [unq_floor,mnout]=calcmeanout(vlscomb)
floorvls=floor(vlscomb(:,1))
unq_floor=unique(floorvls);
for ii=1:length(unq_floor)
    ind=find(floorvls==unq_floor(ii));
    mnout(ii,2)=mean(vlscomb(ind,2));
    mnout(ii,1)=unq_floor(ii)+.5
    
end

function [unq_floor,synout]=calcsynout(vlscomb)
floorvls=floor(vlscomb{1}(:,1))
unq_floor=unique(floorvls);
for ii=1:length(unq_floor)
    for ntind=1:length(vlscomb)
        crvls=vlscomb{ntind}
        ind=find(floor(crvls(:,1))==unq_floor(ii));
        synout(ii,ntind)=length(ind);
    end
end
  