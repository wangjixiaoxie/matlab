

%inactiv_fig2script.m
%main panel is made with call to plothists, and call to plotline.

%raw data panel is made with calls to inactivplotoneday.m
%rewritten 10.11.09 in order 

%modified to plot pre and post acsf histograms separately.
%redesigning this to show first black, blue, grey- baseline run, with
%differnt widths.
%then pitch shift acsf am.
%Then inactiv run

%plothistinactivtw3,
%new figure with basic baseline effect.
function [ax]=plothistinactivtw6(avls,ps)


col{1}=[0 0 0]
col{2}=[0.4 0.4 1]
col{3}=[0.5 0.5 0.5]
ax=ps.ax;
nts=ps.nts;
muind=ps.muind;
initmean=ps.initmean;
plotepochs=ps.plotepochs
sfact(1)=3000
numhist=ps.numhist;
% sfact(2)=1000;

% ax(1)=subplot(511)
% ax(2)=subplot(512)
% ax(3)=subplot(513)
% ax(4)=subplot(514)
% ax(5)=subplot(515)
% ax(6)=subplot(236)

%total traces
%9 traces total.
%first subplot, just acsf
%second subplot, just acsf
%third subplot, three traces
%fourth subplot, one trace
%fifth subplot, three traces


% acmuvec=[1 2 3  1  1 2 3  1   1 2 3 ];
ntvec=[ones(1,3)*nts(1)]
% muvec=[1 1 1  2  2 2 2  3  3 3 3 ]
colvec=[1 2 3]

% for ii=1:length(ps.plotepochs)
% %     if(ii<9)
%         initmean(ii)=avls.acmean(1,avls.mulist(:,muind(1)));
% %     else
% %         initmean(ii)=avls.acmean(2,avls.mulist(:,muind(1)));
% %     end
% end
for ii=1:numhist
    hs(ii).plothist=1;
    hs(ii).ax=ax;
    hs(ii).wid=[2];
    hs(ii).verline=initmean;
    hs(ii).col=col{colvec(ii)};
    hs(ii).orient=['norm'];
    hs(ii).acvec=plotepochs(ii);
    hs(ii).ntvec=ntvec(ii);
    hs(ii).muvec=avls.mulist(muind(ii));
%     hs(ii).muvec=avls.mulist(muind(muvec(ii)));
    %first 6 values are ectarget note
    hs(ii).ntvl=nts(1);
    hs(ii).horline=[];
    hs(ii).verline=hs(ii).verline/sfact(1);
        %second 6 values are control note.
    hs(ii).lw=3;
    hs(ii).ls='-'
end

%assign hstvals for plot
% loadsumdata(bs,2);
hs=gethstcomb2(hs,avls);
for ii=1:numhist
    if(ii<12)
        hs(ii).edges=hs(ii).edges./sfact(1)
    else
        hs(ii).edges=hs(ii).edges./sfact(2)
    end
end
plothists2(hs);

%also create a linestruct which plot 18 lines, one for each acsf/mu run.
%each line needs a color
% lns(ii).col
% ls(ii).ax
% ls(ii).ypts
% ls(ii).xpts
% ls(ii).lw
% ls(ii).ls

%call plotlines2
%count up axes
%first set up the initial 12 structures.

% colvec=[1 2 1 1 2 1 1 2 1 1 2 1 3 1 3 1]
% axvec= [1 1 1 3 3 3 4 4 4 6 6 6 2 2 5 5]
% vlsvec=[4 5 6 9 10 12 4 5 6 9 10 12 4 9 4 9]
% ntvec= [1 1 1 1 1 1 2 2 2 2 2 2 1 1 2 2]
% lns=[];
% for ii=1:16
%     nt=ntvec(ii);
%     mus=avls.mulist(muind(muvec(ii)));
%     if(acmuvec(ii))
%         mn=avls.acmean(nt,mus);
%         stdv=avls.acstdv(nt,mus);
%     else
%         mn=avls.mumean(nt,mus);
%         stdv=avls.mustdv(nt,mus);
%     end
%     if(mod(ii,2)==1)
%         xvl=-0.1;
%    
%     else
%         xvl=-0.1;
%     end
%     marklocwid=.05
%     lns(ii).col=col{colvec(ii)};
%     lns(ii).ax=ax(ceil(ii/2));
%     lns(ii).xpts=[xvl xvl]
%     lns(ii).ypts=[mn mn]
%     if(ii<12)
%         lns(ii).ypts=lns(ii).ypts./sfact(1);
%     else
%         lns(ii).ypts=lns(ii).ypts./sfact(2);
%     end
%     lns(ii).lw=1;
%     lns(ii).ls='-'
%     lns(ii).plotmarker=1;
%     
%     lns(ii).markloc=[xvl-marklocwid xvl mean(lns(ii).ypts) mean(lns(ii).ypts)];
%     lns(ii).marksize=4;
%    
% end

%now 

% plotlines2(lns);

% tsts=[2 3 4 ]
% for ii=1:length(tsts)
%     ind=tsts(ii);
%     axes(ax(ind));
%     box off;
%     axis off;
% end
% 
% linkaxes(ax(1:4));
% linkaxes(ax(5:8));
% % axes(ax(1));
% % axis([0 0.8 2.15 2.4])
% % axes(ax(4));
% % axis([0 0.8 2 2.6])
% 
% 
% %two extra bits, put in wn red line at 2.24
% axes(ax(2));
% plot([0 0.5], [2.245 2.245],'r--');
% %put in arrow for arrowhead.
% axes(ax(3))
% plot([avls.acmean(1,10)/3000 avls.mumean(1,10)/3000],[.75 .75],'k-','Linewidth',2)


%set up the final panel
% exps.plotextraline=0;
% exps.indtoplot=10;
% exps.ntind=1;
% exps.col={[0 0 0]  [0 0 0]  [0.4000 0.4000 1]  [0.4000 0.4000 1]}
% exps.runtype='ex'
% exps.birdind=2
% exps.ax=subplot(1,8,8);

% inactivplotcomparsingle(avls,exps,bs)



% linkaxes([hs(1:3).ax, lns.ax],'y'); 


% %plot scale bar
% axes(lns.ax);
% plot([-.2 -.2+2/24], [6600 6600],'Linewidth',5,'color','k');
%  axis off;
% axes(hs(4).ax)
% axis off;
% axes(exps.ax)
% axis([0.49 0.51 6500 7100])
