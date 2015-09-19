
%option 1.
%edited so that it plots wn, nighttime, raw points, and morning evening


function [os_ph,os_st,plotvls_st,offz]=plot_consolid1TEST(phsumbs,sumbs,bs,sts)
figure

% figstoplot=[1 2 3 4 5  ]
% figstoplot=[5 6 7 ];
% figstoplot=[1 2 7 12]
% figstoplot=[1 2 3 4 7 8 12 13 ]
% figstoplot=[1 2 3 4 7 8 12 13 ];
%EVERYTHING BUT REVERSE 
figstoplot=[ 4 5 6 7 ];
% figstoplot=1;
% figstoplot=[4]
% settings for first figure (time course)


  avlspath='/oriole/bk48w74/datasum'
  avlsmat='datsum.mat'
    
  cmd=['cd ' avlspath];
  eval(cmd);
  cmd=['load ' avlsmat];
  eval(cmd);

ps(1).bsnum=[2 1]
ps(1).runstoplot{1}=[1 4 9 10 11 12 13 14 15 16 17]
ps(1).runstoplot{2}=[2 3 4 5 6 7 13  14 15 16 17]

ps(2).bsnum=[2 ]
ps(2).runstoplot{1}=[4 13]
% ps(2).runstoplot{2}=[3 14]
ps(2).totalplot=2;
ps(2).shiftdrxn=[1 0]
ps(2).arrowht=[0.4 0.45]
ps(2).distht=0.35
ps(2).shiftnum=1;

crbs=2;
    
    cmd=['cd ' bs(crbs).path ];
    eval(cmd);
    cmd=['load ' bs(crbs).matfilename];
    eval(cmd);
clear ps;
%Example asymptote plot.
% figure
%   [sumdyn,sumdynasymp,sumdynrev]=selectdynamicruns3(sumbs,0,0,0,0)


if (ismember(3,figstoplot))
    crbs=2;
    
    cmd=['cd ' bs(crbs).path ];
    eval(cmd);
    cmd=['load ' bs(crbs).matfilename];
    eval(cmd);
    
%     ps.STIM=0;
        ps.plotarrow=1;
%         ps.runstoplot=[5 7 8 9 10 11 12 13 14 15]
        ps.runstoplot=[1:11 ]
         ps.wncol=[.98 .8 .9]
%         crax=subplot(5,5,1:4)
        crax=subplot(3,2,1:2)        
        ps.ax=crax;
        ps.muht=2.6;
        ps.plotsep=1;
        ps.adjx=1;
        
        ps.acindvl=[14:19]
        ps.muindvl=[6:9]
        ps.FILL=0;
         ps.plotedgeybnds=[2 3]
         ps.startind=13;
         ps.endind=15;
        ps.plotedgexbnds=[0:15]
        ps.plotwn=1;
        ps.plotz=0
        ps.wnfloor=2.2
        ps.plotwnline=1;
        ps.plotmumarker=0;
        ps.plotmuline=0;
        ps.plotacline=0;
        ps.plotrawmu=1;
        ps.plotrawac=0;
        ps.plotallac=1;
       
        ps.plotallacpts=1;
        ps.plotinitendline=0;
        ps.flip=0;
        ps.mumarksize=1.5
%         axtmcourse(ii)=subplot(1,bsln,ii);
        plot_tmcourse5MOD(phsumbs(crbs),ps,avls);
%         axis square
        axis([ -2 9 2.25 2.7])
%         ps.plotedge=1;
%         ps.plotedgeybnds=[2 3]
%         ps.plotedgexbnds=[0:15]
end

if(ismember(4,figstoplot))
crax=subplot(3,4,1:2)
    ps.ax= crax
          ps.mulistind= [1 2 3 4 5 6 7 8 9 10 11]
               ps.STIM= 0
          ps.plotarrow= 0
            ps.plotave= 1
       ps.plotmumarker= 0
            ps.acindvl= [10:18]
            ps.muindvl= [10 11 13 15:18]
           ps.startind= 1
              ps.plotz= 0
               ps.flip= 1
               ps.startind=phsumbs(2).aclist(11,1);
               ps.endind=15;
               ps.FILL= 0
            ps.wnfloor= 0
         ps.plotwnline= 1
              ps.wncol= 'r'
            ps.plotsep= 1
            ps.plotdiff=1;
    ps.plotinitendline= 0
         ps.plotmuline= 1
          ps.plotrawmu= 0
          ps.plotallac= 1
          ps.netdiff=1;
          ps.plotallacpts=0;
          ps.plotrawac= 0
         ps.runstoplot= [10:18]
              ps.muind= [1 2 3 4 5 6 7]
         ps.plotacline= 1
       
         [os_ph]=plot_tmcourse5MOD2(phsumbs(2),ps,avls)
        tst=1;
end

if(ismember(5,figstoplot))
ptsps.ax=subplot(3,4,3:4)
   
    ptsps.axbnds=[-1 9 2.25 2.55]
    ptsps.mulistind=1;
    ptsps.plotind=[1:7 9:12]
    ptsps.initind=sumbs(2).mulist(2);
    ptsps.STIM=1;
    ptsps.plotarrow=0;
    ptsps.netdiff=1;
    ptsps.plotdiff=1;
    ptsps.startind=1;
    [os_st]=plot_tmcourse2rev(sumbs(2),ptsps)
   axis([-1 9 2.25 2.55])
end
%plot fill
if(ismember(6,figstoplot))
    ps.ax=subplot(3,4,5:6);
    mkfillnetdiff(os_ph,ps);
    axis([-2 11 -80 80])
    axis
    ps.ax=subplot(3,4,7:8);
    ps.flip=0;
    mkfillnetdiff(os_st,ps);
    axis([-2 11 -80 80])
end

if(ismember(7,figstoplot))
    axin(1)=subplot(3,4,9:10);
    axin(2)=subplot(3,4,11:12);
    [plotvls_st,offz]=plothistoryfigscript6(sumbs,phsumbs,axin);
end
function []=mkfillnetdiff(os,ps)
    
fillx=[makerow(os.nettms) makerow(os.nettms(end:-1:1))];
    filly=[os.netdiff zeros(1,length(os.netdiff))];
col=[0.4 0.4 0.4]
    
hold on;
if(~ps.flip)
plot(os.nettms,os.netdiff,'ko','Linewidth',1,'MarkerSize',6,'MarkerFaceColor','k');
fill(fillx,filly,col,'EdgeColor','k');
else
    plot(os.nettms,-os.netdiff,'ko','Linewidth',1,'MarkerSize',6,'MarkerFaceColor','k');
    fill(fillx,-filly,col,'EdgeColor','k');
end