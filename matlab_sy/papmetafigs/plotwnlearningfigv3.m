function [vlscomb]=plotwnlearningfigv3(vlscomb)
%first arm up vls file

figstoplot=[ 1 2 4];
if(~exist('vlscomb'))

%aug 18 and aug19
pth{1}='/oriole2/bk57w35evren/bk57w35/run2/stimoffaug18'
cvl{1}='batch.rand.keep'



pth{2}='/oriole2/bk57w35evren/bk57w35/amponaug19'
cvl{2}='batchcomb'

%aug22
pth{3}='/oriole2/bk57w35evren/bk57w35/amponaug19'
cvl{3}='batch22.randb.rand'

pth{4}='/oriole2/bk57w35evren/bk57w35/amponaug19/23files'
cvl{4}='batch23.rand'

pth{5}='/oriole2/bk57w35evren/bk57w35/amponaug19/23_24_25files'
cvl{5}='batch24.rand'
pth{6}='/oriole2/bk57w35evren/bk57w35/amponaug19/23_24_25files'
cvl{6}='batch25.rand'
pth{7}='/oriole2/bk57w35evren/bk57w35/amponaug19/23_24_25files'
cvl{7}='batch26.rand'
pth{8}='/oriole2/bk57w35evren/bk57w35/amponaug19/27files'
cvl{8}='batch27.rand'

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
  tbinshft=.084
  NFFT=512;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
  fbins=[6000 8000];
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
    
    [vals,trigs]=triglabel(bt,'a',1,1,0,0);

end
end

if(ismember(1,figstoplot))
   subplot(3,4,9:12)
   [unq_floor,mnout]=calcmeanout(vlscomb);
   startday=3;
   init_tm=floor(mnout(startday,1));
   plot(vlscomb(:,1)-init_tm,vlscomb(:,2)/3000,'ko','MarkerSize',0.8,'MarkerFaceColor','k');
   hold on;
   indht=find(vlscomb(:,4)>0&vlscomb(:,4)<105&(vlscomb(:,1)-init_tm>0))
   
   plot(vlscomb(indht,1)-init_tm,vlscomb(indht,2)/3000,'ro','MarkerSize',0.8,'MarkerFaceColor','r');
   plot(mnout(:,1)-init_tm,mnout(:,2)/3000,'ko','MarkerSize',4)
    plot(mnout(:,1)-init_tm,mnout(:,2)/3000,'k','Linewidth',1)
    daspect([30 1 1])
end

%plot three histograms - init histogram over last baseline day.
%day 3 all
%day 3 hits
if(ismember(2,figstoplot))
    subplot(3,4,7:8);
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
    plot(mnbas,0.3,'Marker','v','Color','k');
    plot(meanht,0.3,'Marker','v','Color','b');
end  
  
if(ismember(3,figstoplot))
    subplot(3,4,7:8);
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
     
    exsong(1).path='/oriole2/bk57w35evren/bk57w35/amponaug19'
%     exsong(1).cvl=cvl{2}
    exsong(1).fn='bk57w35_200809_083338.2512.cbin';
    exsong(1).bnds=[4.2 6.2]
    exsong(1).ax=subplot(3,4,3:4)
    ps.detect=0;
    plotcbin(exsong);
    
    exsong(1).path='/oriole2/bk57w35evren/bk57w35/run2/stimoffaug18'
%     exsong(1).cvl=cvl{2}
    exsong(1).fn='bk57w35_180809_073117.7369.cbin';
    exsong(1).bnds=[3.2 5.2]
    exsong(1).ax=subplot(3,4,1:2)
    ps.detect=0;
    plotcbin(exsong);
    
    
    exsong(1).path='/oriole2/bk57w35evren/bk57w35/run2/stimoffaug18'
%     exsong(1).cvl=cvl{2}
    exsong(1).fn='bk57w35_180809_073117.7369.cbin';
    exsong(1).bnds=[3.55 3.7]
    exsong(1).ax=subplot(3,4,5:6)
    ps.detect=0;
    plotcbin(exsong);
    hold on;
    plot([.0562 .0562],[0 10000],'r--');
    hold on;
     plot([.0722 .0722],[0 10000],'r--');
    
%     exsong(1).stim_path='/oriole/bk48w74/stimexample'
%     exsong(1).stim_mat='stimexample.mat'
end
    

function [unq_floor,mnout]=calcmeanout(vlscomb)
floorvls=floor(vlscomb(:,1))
unq_floor=unique(floorvls);
for ii=1:length(unq_floor)
    ind=find(floorvls==unq_floor(ii));
    mnout(ii,2)=mean(vlscomb(ind,2));
    mnout(ii,1)=unq_floor(ii)+.5
    
end


  