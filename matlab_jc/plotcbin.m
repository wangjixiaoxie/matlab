%example inputs.

% exsong.path=avls.pvls{16}
% exsong.fn='bk63w42_051208_1706.122.cbin';
% exsong.bnds=[5.7 7.3]
% exsong.boxtoplot=[4:6]
% exsong.scalfreqs=[6000 8000]
% 
% plotboxes(1).ntnm=1;
% plotboxes(2).ntnm=2;
% plotboxes(3).ntnm=3;
% plotboxes(4).ntnm=1;
% plotboxes(5).ntnm=2;
% plotboxes(6).ntnm=3;
% 




function [datout,fs,ax]=plotcbin(exsong, ps)
ps.NORMFREQS=1;
axes(exsong.ax);
ax=exsong.ax
%     axes(ax);

path=exsong.path;
bnds=exsong.bnds;
fn=exsong.fn;
PLOT_TRIG=0
strcmd=['cd ' path]
eval(strcmd);
[plainsong,fs]=ReadCbinFile(fn);
% if(isfield(ps,'plot_trigs'))
%     if(ps.plot_trigs)
%         rd=readrecf(fn);
%         PLOT_TRIG=1;
%     end
% end
%here for wav, =wavread('');

datout=plainsong(floor(bnds(1)*fs):floor(bnds(2)*fs),:);
plainsong=floor(plainsong(bnds(1)*fs:bnds(2)*fs));
[sm,sp,t,f]=evsmooth(plainsong,fs,500);

if(exist('normfreqs'))
    spt{1}=log(abs(sp))
    ps.freqlimbound=[6000 8000]
    ps.freqspacing=62.5
    [colbnds]=calcbndsavn(spt,ps)
    [colbnds]=1.3*colbnds;
end


 imagesc(t+.008,f,log(abs(sp)));syn;ylim([0,1e4]);
if (exist('normfreqs'))
 set(gca,'CLim',colbnds);
end
 hold on;
 if(exist('plotlns'))
    xvec=[plotlns-bnds(1);plotlns-bnds(1)];
    yvec=[zeros(length(plotlns),1)';10000*ones(length(plotlns),1)'];
    plot(xvec,yvec,'w')
    
 end
 if(PLOT_TRIG)
     hold on;
     plot(rd.ttimes/1000-bnds(1),9000*ones(length(rd.ttimes),1),'v','MarkerFaceColor','c','MarkerEdgeColor','c','MarkerSize',5);
 end
 
 if (exist('plotboxes'))
%      1:length(plotboxes)
for ii=1:length(plotboxes)
        plotboxes(ii).fn=fn;
        
        plotboxes(ii).bnds=bnds;
        [vls]=getbnds(plotboxes(ii),bnds);
        plotbox(vls,ax)
    end
    end
 colormap('hot')
maxnorm=max(abs(plainsong));
% maxnorm=maxnorm/25;

function [ax]=plotbox(vls,ax)
    col='c';
    xbnds=vls(1:2)
    ybnds=vls(3:4);
%     axes(ax);
    %plot 4 lines;
    plot([vls(1) vls(2)],[vls(3) vls(3)],'Color',col)
    plot([vls(1) vls(2)],[vls(4) vls(4)],'Color',col)
    plot([vls(1) vls(1)],[vls(3) vls(4)],'Color',col)
    plot([vls(2) vls(2)],[vls(3) vls(4)],'Color',col)
function [xyvc]=getbnds(ps,bnds)
    ntnm=ps.ntnm
    sfact=ps.sfact;
    ofnm=ps.lbl;
    yvls=ps.fbins/sfact;
    tshft=ps.tshft
    NFFT=ps.NFFT
    fnn=[ps.fn,'.not.mat'];
    
    
    load(fnn);
    ofst_tm=onsets(ofnm)/1000-bnds(1);
    xvls(1)=ofst_tm+tshft
    xvls(2)=ofst_tm+tshft+NFFT/44100;
    xyvc=[xvls yvls]

    
% 
% nrmsong=32768*(1/(1.01*maxnorm))*plainsong;
% 
% %writeaiff_file
% aif_file='example.aiff'
% aiffwrite(aif_file,nrmsong,32000,16);
% 
%  
%  %  figure
% %  if (labels)
% %      load([fnvl{ii} '.not.mat'])
% %      zoomind=find(onsets/1000>bnds{ii}(1)&onsets/1000<bnds{ii}(2));
% %      for jj=1:length(zoomind)
% %          text(onsets(zoomind(jj))/1000-bnds{ii}(1), 200, labels(zoomind(jj)),'Color', 'w','Fontsize',16);
% % % ,'Color','w');
% %      end
%  end
%  box off;
%  axis on;   
% 
% outwav=wavwrite(nrmsong,32000, 16,'example.wav')
