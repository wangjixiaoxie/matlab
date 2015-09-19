%REVISED 8.4.09

%create baseline contour.
%modified 7.27.09 in order to make multiple plots;
% stimx{1}=[.055 .07]
% stimx{2}=[.04 .08]
% stimx{3}=[.065 .08]
% muinds=[1:24]
PLOT_TMCOURSE=1;
HS_TO_PLOT=[1:47]
% HS_TO_PLOT=[1:43]
clear mnvlsfb mnvlsct stdvlsfb stdvlsct tmvec htvlsct htvlsfb axcnt

% cd /oriole6/bk59w37redo/wnon710
% 
% load pitchdata.mat
% pitchdata_bas=contours;
% meanbas=mean(pitchdata_bas,2);
% stdbas=std(pitchdata_bas,0,2);
figure
%  muinds=[4 17:19]`
PLOTCONTOURS=1;
figsperplot=8;
%for bk57w35
TIMEBNDS=[.07 .12]

% %for r34w20
% TIMEBNDS=[.02 .07]
% 
%BELOW is bk48w74 long pulse
% muinds=[1:7 13:15 17:19 23 24 25 28 29 31 33 34 35]

muinds=[5:27]


% muinds=[1:40 42 44 45]
ln=length(muinds);

%THIS LOOP MAKES CONTOURS
if(PLOTCONTOURS)
    for ii=1:length(muinds)
    ii
    crind=muinds(ii);
    pathvl=avls.pvls{crind}
     btvl=avls.cvl{crind};
    
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]
    else
        cmd=['cd '  pathvl]
    end
    eval(cmd);
    cmd=['load ' btvl '.mat']
    eval(cmd);
    
    modvl=mod(ii,figsperplot)
    
%     if(modvl==0)
%        figure; 
%     end
%     plotnum=modvl+1;
%     axcnt(ii)=subplot(figsperplot,1,plotnum);
%     if(~isempty('crfbind'))
%         mnfb=mean(contours(:,crfbind),2);
%     end
%     mnct=mean(contours(:,crctind),2);
%     if(ii==1)
%         basct=mnct
%     end
%     if(~isempty('crfbind'))
%      
%         stdfb=std(contours(:,crfbind),0,2);
%         stefb=stdfb./sqrt((length(crfbind)));
% %     end
%         stdct=std(contours(:,crctind),0,2);
%     stect=stdct./sqrt(length(crctind));
%     ind=find(pitchtms>TIMEBNDS(1)&pitchtms<TIMEBNDS(2))
%     if(~isempty('crfbind'))
%         plot(pitchtms(ind), mnfb(ind)-stefb(ind),'r');
%         hold on;
%         plot(pitchtms(ind),mnfb(ind),'r','Linewidth',2)
%         plot(pitchtms(ind), mnfb(ind)+stefb(ind),'r');
%     end
        
        plot(crind, vals(crctind,2),'k.');
        hold on;
        plot(crind, vals(crfbind,2),'r.');
    end    
end
% figure
% %THIS LOOP MAKES HISTOGRAMS
% for ii=1:length(HS_TO_PLOT)
%     crind=muinds(HS_TO_PLOT(ii));
%     pathvl=avls.pvls{crind}
%     btvl=avls.cvl{crind};
%     if(isfield(avls,'baspath'))
%        cmd=['cd ' avls.baspath pathvl]
%        eval(cmd);
%     else
%         cmd=['cd ' pathvl]
%         eval(cmd);
%     end
%     cmd=['load ' btvl '.mat']
%     eval(cmd);
% %     load extradata.mat
%       modvl=mod(ii,figsperplot)
%      if(modvl==0)
%        figure; 
%   
%     end
%     plotnum=modvl+1;
%     
%     axhist(ii)=subplot(figsperplot,1,plotnum);
%    stairs(HST_EDGES/3,hsctnrm,'k')
%    
%    hold on;
%    if(ii==1)
%        mnbas=ctmean/3;
%        stdbas=stdct/3;
%    end
%    
% %    plot([2294 2294],[0 1],'c--','Linewidth',2)
%    plot([mnbas+stdbas mnbas+stdbas], [0 1], 'c--')
%    plot([mnbas-stdbas mnbas-stdbas], [0 1], 'c--')
%    plot([mnbas mnbas],[0 1],'c--','Linewidth',2)
%    tmvec(ii,:)=tms;
%    mnvlsct(ii)=ctmean; 
%    stdvlsct(ii)=stdct;
%    text(2125,0.5,['catch=' num2str(length(crctind))],'Color','k');
%    plot([(mnvlsct(ii)-stdvlsct(ii))/3 (mnvlsct(ii)+stdvlsct(ii))/3],[0.4 0.4],'k')
%    
%    if(~isempty(crfbind))
%         mnvlsfb(ii)=fbmean;
%         stdvlsfb(ii)=stdfb;
%         plot([(mnvlsfb(ii)-stdvlsfb(ii))/3 (mnvlsfb(ii)+stdvlsfb(ii))/3],[0.45 0.45],'r')
%         text(2125,0.2,['fb=' num2str(length(crfbind))],'Color','r');
%         stairs(HST_EDGES/3,hsfbnrm,'r')
%         box off;
%    end
% end
%     
%   
%     
%       
% 
% 
% if(PLOT_TMCOURSE)
%     mnvlsfb=[mnvlsfb]
%     mnvlsct=mnvlsct
%     htvlsfb=stdvlsfb
%     htvlsct=stdvlsct
%     plot_tmcourse_general(tmvec,mnvlsfb,mnvlsct,htvlsfb,htvlsct);
%     
% end
% 
% 
%     
%     