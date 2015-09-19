%plotstimfig1.
clear catchvl

%%%this is spectrogram component.

exsong(1).path=[avls.baspath avls.pvls{2}]
exsong(1).fn='bk48w74_210809_105358.168.cbin';
exsong(1).bnds=[4 6.7]
exsong(1).fvst=fvstex{1};

exsong(2).path=[avls.baspath avls.pvls{4}]
exsong(2).fn='bk48w74_230809_180859.15.cbin';
exsong(2).bnds=[7.3 10]
exsong(2).fvst=fvstex{2}
clear ax

figure;
ax1(1)=subplot(6,1,2:6)
exsong(1).ax=ax1(1);
ax1(2)=subplot(6,1,1);

figure
ax2(1)=subplot(6,1,2:6);
exsong(2).ax=ax2(1);
ax2(2)=subplot(6,1,1);
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
for ii=1:2

    [dat,fs]=plotcbin(exsong(ii));
    hold on;
    plot([0 2.5],[6600 6600],'c--')
end
% 
%     tms=0:1/fs:length(dat(:,2))/fs
%     tms2=tms(1:end-1);
%     ax2(2)=subplot(6,1,1)
%     plot(tms2,dat(:,2),'k','Linewidth',2)
axes(ax1(2));
%this is for first one.    
fvst=fvstex{1}
    onstime(1)=fvst(2).ons(fvst(2).ind)
onstime(2)=fvst(3).ons(fvst(3).ind)
onstime(3)=fvst(4).ons(fvst(4).ind)
% onstime(3)=fvst(9).ons(fvst(9).ind)
stimofftime(1)=fvst(2).STIMTIME;
catchvl(1)=fvst(2).STIMCATCH
stimofftime(2)=fvst(3).STIMTIME;
catchvl(2)=fvst(3).STIMCATCH
stimofftime(3)=fvst(4).STIMTIME;
catchvl(3)=fvst(3).STIMCATCH;
dtime(1)=onstime(1)/1000+stimofftime(1)/1000-exsong(1).bnds(1);
dtime(2)=onstime(2)/1000+stimofftime(2)/1000-exsong(1).bnds(1);
dtime(3)=onstime(3)/1000+stimofftime(3)/1000-exsong(1).bnds(1);

hold on;
for ii=1:3
    plot([dtime(ii) dtime(ii)],[0 20000],'r','Linewidth',2)
end

for ii=1:3
    if catchvl(ii)<1
    plot(dtime(ii)+.065+stimtms,stimblock,'k','Linewidth',2)
    end
end
    
%this is for second one.
axes(ax2(2));
fvst=fvstex{2}
%these happen to be 19 and 20.
onstime(1)=fvst(40).ons(fvst(40).ind)
onstime(2)=fvst(41).ons(fvst(41).ind)
onstime(3)=fvst(42).ons(fvst(42).ind)
% onstime(3)=fvst(9).ons(fvst(9).ind)
stimofftime(1)=fvst(40).STIMTIME;
catchvl(1)=fvst(40).STIMCATCH
stimofftime(2)=fvst(41).STIMTIME;
catchvl(2)=fvst(41).STIMCATCH;
stimofftime(3)=fvst(42).STIMTIME;
catchvl(3)=fvst(42).STIMCATCH
dtime(1)=onstime(1)/1000+stimofftime(1)/1000-exsong(2).bnds(1);
dtime(2)=onstime(2)/1000+stimofftime(2)/1000-exsong(2).bnds(1);
dtime(3)=onstime(3)/1000+stimofftime(3)/1000-exsong(2).bnds(1);

hold on;
for ii=1:3
    plot([dtime(ii) dtime(ii)],[0 20000],'r','Linewidth',2)
end

for ii=1:3
    if catchvl(ii)<1
    plot((dtime(ii)+.065+stimtms),stimblock,'k','Linewidth',2)
    end
end


title('plotstimexample.m')

%%%%this is the part which shows example time course of pitch.
%%%plot two panels, 
%%%one panel for initial downshift, other panel for upshift.
%%%PLOTPANEL2.m



% 
% 
% linkaxes(ax2,'x')
% linkaxes(ax1,'x')


%contour figure
% 
% PLOTCONTOURS=1;
% muinds=[45 55]
% if(PLOTCONTOURS)
%     for ii=1:length(muinds)
%     
%     crind=muinds(ii);
%     pathvl=avls.pvls{crind}
%      btvl=avls.cvl{crind};
%     
%     if(isfield(avls,'baspath'))
%         cmd=['cd ' avls.baspath pathvl]
%     else
%         cmd=['cd '  pathvl]
%     end
%     eval(cmd);
%     cmd=['load ' btvl '.mat']
% 
%     eval(cmd);
%     
%     modvl=mod(ii,figsperplot)
%     
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
%     end
%         stdct=std(contours(:,crctind),0,2);
%     stect=stdct./sqrt(length(crctind));
%     ind=find(pitchtms>TIMEBNDS(1)&pitchtms<TIMEBNDS(2))
%     if(~isempty('crfbind'))
%         plot(pitchtms(ind), mnfb(ind)-stefb(ind),'r');
%         hold on;
%         plot(pitchtms(ind),mnfb(ind),'r','Linewidth',2)
%         plot(pitchtms(ind), mnfb(ind)+stefb(ind),'r');
%     end
%         
%         plot(pitchtms(ind), mnct(ind)+stect(ind),'k');
%         plot(pitchtms(ind), mnct(ind)-stect(ind),'k');
%         plot(pitchtms(ind),mnct(ind),'k','Linewidth',2)
%         if(ii>1)
%             plot(pitchtms,basct,'c--')
%         end
% 
%     end
%     linkaxes(axcnt(1:length(muinds)));
%     box off;
%     
% end









































































































































































































































































































































































































































































































































