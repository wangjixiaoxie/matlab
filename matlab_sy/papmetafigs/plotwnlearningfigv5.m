function [vlscomb,sumdata]=plotwnlearningfigv4(vlscomb)
%first arm up vls file

figstoplot=[ 1:5];

threshdays=[0 1 1 2 2 8]
threshvals=[2.44 2.44 2.47 2.47 2.52 2.52]
if(~exist('vlscomb'))

%oct 29 and oct30
pth{1}='/oriole/bk48w74/BASELINE4/baseline1028'
cvl{1}='batch29.rand.keep'



pth{2}='/oriole/bk48w74/BASELINE4/temptest/files30'
cvl{2}='batchcomb'

%aug22
pth{3}='/oriole/bk48w74/BASELINE4/1031'
cvl{3}='batchcomb.rand.rand'

pth{4}='/oriole/bk48w74/BASELINE4/temptest/files01'
cvl{4}='batch.rand.rand.keep'

pth{5}='/oriole/bk48w74/BASELINE4/temptest/files02'
cvl{5}='batchcomb'
pth{6}='/oriole/bk48w74/BASELINE4/temptest/files03'
cvl{6}='batch.rand.keep'

pth{7}='/oriole/bk48w74/BASELINE4/temptest/files04'
cvl{7}='batch.rand.keep'

pth{8}='/oriole/bk48w74/BASELINE4/temptest/files05'
cvl{8}='batch.rand.keep'
pth{9}='/oriole/bk48w74/BASELINE4/temptest/files06'
cvl{9}='batch.rand.keep'

pth{10}='/oriole/bk48w74/BASELINE4/temptest/files07'
cvl{10}='batch.rand.keep'


%notes day 1 thresh varied from 2.42 to 2.46, plot as 2.44


% pth{9}='/oriole2/bk57w35evren/bk57w35/amponaug19'
% cvl{9}='batch28.rand'
vlscomb=[];
SOUNDSTR='wav'
STIMSTR='rig'
STIMBNDS=[200 -200]
for ii=1:length(pth)    
  cmd=['cd ' pth{ii}];
  eval(cmd);
  bt=cvl{ii};
  tbinshft=.074
  NFFT=512;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
  fbins=[6000 8300];
  save BINS_B NFFT fbins tbinshft
  % frequency analysis just for 'b'
  load BINS_B
  NT='a';PRENT='';PSTNT='';
    fv=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
     fvst=findwnote9(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0',0,STIMBNDS, SOUNDSTR,STIMSTR);
%     fvd=findwnote4(bt2,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
    vls=getvals_sec(fv,1,'TRIG');
for ii=1:length(fvst)
       if(~isempty(fvst(ii).SOUNDTIME))
           vls(ii,4)=fvst(ii).SOUNDTIME;
       else
           vls(ii,4)=0
       end
    end
    
    vlscomb=[vlscomb;vls]
    sumdata=[];
    
    [vals,trigs]=triglabel(bt,'a',1,1,0,0);

end
end

if(ismember(1,figstoplot))
   figure
   sumdata=[];
   subplot(3,4,5:8)
   [unq_floor,mnout,stdout]=calcmeanout(vlscomb);
   startday=3;
   init_tm=floor(mnout(startday,1))+.5;
   plot(vlscomb(:,1)-init_tm,vlscomb(:,2)/3000,'o','MarkerSize',1,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor','none');
   hold on;
   indht=find(vlscomb(:,4)>0&vlscomb(:,4)<105&(vlscomb(:,1)-init_tm>0))
   
   plot(vlscomb(indht,1)-init_tm,vlscomb(indht,2)/3000,'o','MarkerSize',1,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor','none');
   plot(mnout(:,1)-init_tm,mnout(:,2)/3000,'ko','MarkerSize',4)
   plot([mnout(:,1)'-init_tm; mnout(:,1)'-init_tm],[(mnout(:,2)/3000)'-stdout/3000; (mnout(:,2)/3000)'+stdout/3000],'k')
    plot(mnout(:,1)-init_tm,mnout(:,2)/3000,'k','Linewidth',1.5)
    daspect([10 1 1])
    plot(threshdays,threshvals,'r--')
end

%plot three histograms - init histogram over last baseline day.
%day 3 all
%day 3 hits
if(ismember(2,figstoplot))
%     figuremnout(:
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
%     subplot(3,4,7:8);
    edges=[6600:40:7600]
    floorvls=floor(vlscomb(:,1));
    indbas=find(floorvls<=unq_floor(2));
    hstbas=histc(vlscomb(indbas,2),edges);
    mnbas=mean(vlscomb(indbas,2))/3000;
    hstbas=hstbas./length(indbas);
    
    ind_three=find(floorvls==unq_floor(5));
    hst_three=histc(vlscomb(ind_three,2),edges);
    indht_three=intersect(ind_three,indht);
    hstht=histc(vlscomb(indht_three,2),edges);
    meanht=mean(vlscomb(ind_three,2))/3000;
    hstlearn=hst_three./length(ind_three);
    hstht=hstht./length(ind_three);
    edges=edges(1:end-1);
    hstbas=hstbas(1:end-1);
    hstlearn=hstlearn(1:end-1);
    hstht=hstht(1:end-1);
    
    stairs(edges/3000,hstbas,'k')
    hold on;
    stairs(edges/3000,hstlearn,'b')
    stairs(edges/3000,hstht,'r')
    plot(mnbas,0.15,'Marker','v','Color','k');
    plot(meanht,0.15,'Marker','v','Color','b');
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
    
    exsong(1).path='/oriole/bk48w74/BASELINE4/1031'
%     exsong(1).cvl=cvl{2}
    exsong(1).fn='bk48w74_311009_175923.24.cbin';
    exsong(1).bnds=[2.1 3.9]
    exsong(1).ax=subplot(3,4,2:4)
    ps.detect=0;
    plotcbin(exsong,ps);
    
    
    exsong(1).path='/oriole/bk48w74/BASELINE4/temptest/files30'
%     exsong(1).cvl=cvl{2}
    exsong(1).fn='bk48w74_301009_170051.397.cbin';
    exsong(1).bnds=[7.4 7.55]
    exsong(1).ax=subplot(3,4,1)
    ps.detect=0;
    plotcbin(exsong,ps);
    hold on;
    plot([.086 .086],[0 10000],'r--');
    hold on;
     plot([.1020 .1020],[0 10000],'r--');
    
%     exsong(1).stim_path='/oriole/bk48w74/stimexample'
%     exsong(1).stim_mat='stimexample.mat'
end
    

function [unq_floor,mnout,stout]=calcmeanout(vlscomb)
floorvls=floor(vlscomb(:,1))
unq_floor=unique(floorvls);
for ii=1:length(unq_floor)
    ind=find(floorvls==unq_floor(ii));
    mnout(ii,2)=mean(vlscomb(ind,2));
    stout(ii)=std(vlscomb(ind,2));
    mnout(ii,1)=unq_floor(ii)+.5
    
end


  