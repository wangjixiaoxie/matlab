%for plotting stim shifts
function []=plotsameday(sumbs);

% bsind=[1 1 1 1 2 2]
bsind(1)=1;
inds{1}=[ 25 28  ]

bsind(2)=1;
inds{2}=[69 70]




bsind(3)=1;
inds{3}=[14 15]

bsind(4)=1;
inds{4}=[35 36]

bsind(5)=1;
inds{5}=[82 83]

bsind(6)=2;
inds{6}=[13 14 12]

% bsind(7)=2;
% inds{7}=[13 14 12 ]


for ii=1:length(bsind)
    ax(ii)=subplot(length(bsind),1,ii);
    plotedge([-1 6],[-1:2]);
    clear tms acmean mumean preinds pretms premumean preacmean pstinds psttms pstacmean pstmumean
    for jj=1:length(inds{ii})
        crbsind=bsind(ii);
        crbs=sumbs(crbsind);
        
        crtmind=inds{ii}(jj);
        crtm=crbs.tmvec(inds{ii}(1),1);
        crflr=floor(crtm);
        tms(jj)=mean(crbs.tmvec(crtmind,1))-crflr;
        mumean(jj)=crbs.muz(crtmind);
        acmean(jj)=crbs.acz(crtmind);
    end
    
    %get preinds from previous day if any
    
    preinds=find(crbs.flrtmvec(crbs.STANRUNS)-(crflr-1)==0);
    preinds=crbs.STANRUNS(preinds);
    if(isfield(crbs,'NOSTIM'))
        stiminds=find(crbs.NOSTIM(preinds)==0)
        preinds=preinds(stiminds);
    end
    
    [s,sortind]=sort(crbs.tmvec(preinds,1));
    preinds=preinds(sortind);
    for jj=1:length(preinds);
        crind=preinds(jj);
        pretms(jj)=mean(crbs.tmvec(crind,1))-crflr
        premumean(jj)=crbs.muz(crind);
        preacmean(jj)=crbs.acz(crind);
        
    end
    
      pstinds=find(crbs.flrtmvec(crbs.STANRUNS)-(crflr+1)==0);
    pstinds=crbs.STANRUNS(pstinds);
    if(isfield(crbs,'NOSTIM'))
        stiminds=find(crbs.NOSTIM(pstinds)==0)
        pstinds=pstinds(stiminds);
    end
    
    [s,sortind]=sort(crbs.tmvec(pstinds,1));
    pstinds=pstinds(sortind);

    for jj=1:length(pstinds);
        crind=pstinds(jj);
        psttms(jj)=mean(crbs.tmvec(crind,1))-crflr
        pstmumean(jj)=crbs.muz(crind);
        pstacmean(jj)=crbs.acz(crind);
        
    end
    %get postinds from previous day if any
    
    
    
    plot(tms,abs(acmean),'k');
    
    hold on;
    plot(tms,abs(acmean),'ko','MarkerSize',4,'MarkerFaceColor','k')
    plot(tms,abs(mumean),'r');
     plot(tms,abs(mumean),'ro','MarkerSize',4,'MarkerFaceColor','r')
     
     if(~isempty(preinds))
     plot(pretms,abs(preacmean),'k');
       plot(pretms,abs(preacmean),'ko','MarkerSize',4,'MarkerFaceColor','k')
    plot(pretms,abs(premumean),'r');
     plot(pretms,abs(premumean),'ro','MarkerSize',4,'MarkerFaceColor','r')
     end
     if(~isempty(pstinds))
        plot(psttms,abs(pstacmean),'k');
       plot(psttms,abs(pstacmean),'ko','MarkerSize',4,'MarkerFaceColor','k')
    plot(psttms,abs(pstmumean),'r');
     plot(psttms,abs(pstmumean),'ro','MarkerSize',4,'MarkerFaceColor','r')
     end
     axis([-1 2 -1 6])
    plot([-1 2],[0 0],'k--','Linewidth',2)
    
    
    box off
%     axis square;
end

 function []=plotedge(ybnds,xvecin)
%            minx=min([makerow(initxvls) makerow(xvls) makerow(enxvls)]);
%            maxx=max([makerow(initxvls) makerow(xvls) makerow(enxvls)]);
%            st_x=floor(minx);
%            max_x=ceil(maxx);
           
           for ii=1:length(xvecin)
               crx=xvecin(ii);
              xbnds(1)=crx-3/24;
              xbnds(2)=crx+6/24;
              xvec=[xbnds xbnds(end:-1:1)]
              yvec=[ybnds(1) ybnds(1) ybnds(2) ybnds(2)]
              fill(xvec,yvec,[0.7 0.7 0.7],'edgecolor','none');
              hold on;
           end