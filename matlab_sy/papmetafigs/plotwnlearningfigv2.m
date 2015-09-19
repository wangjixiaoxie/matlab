function [combvls]=plotwnlearningfigv2()
figure
figstoplot=[1 3  4 5 ]

  avlspath='/oriole6/bk61w42/datasum'
  avlsmat='pathvals1-analdata.mat'
    
  cmd=['cd ' avlspath];
  eval(cmd);
  cmd=['load ' avlsmat];
  eval(cmd);

  %example spectrum
  if(ismember(1,figstoplot))
     
    exsong(1).path=[ avls.pvls{4}]
    exsong(1).cvl=avls.cvl{4}
    exsong(1).fn='bk63w42_011208_0757.3471.cbin';
    exsong(1).bnds=[5 7.1]
    exsong(1).ax=subplot(421)
    ps.detect=0;
    plotexspect(exsong,ps);
%     exsong(1).stim_path='/oriole/bk48w74/stimexample'
%     exsong(1).stim_mat='stimexample.mat'
end

%insert histogram plot showing histograms stable over time

%then plot showing that you can add whitenoise to histograms

%spectrogram with white noise example
if(ismember(3,figstoplot))
     
    exsong(1).path=[ avls.pvls{16}]
    exsong(1).cvl=avls.cvl{16}
    exsong(1).fn='bk63w42_051208_0754.2248.cbin';
    exsong(1).bnds=[3 5.1]
    exsong(1).ax=subplot(425)
    ps.detect=1;
    plotexspect(exsong,ps);
%     exsong(1).stim_path='/oriole/bk48w74/stimexample'
%     exsong(1).stim_mat='stimexample.mat'
end

%then change in histograms  before and after white noise.
if(ismember(4,figstoplot))
    pvlindpre=13;
    pvlindpst=42;
    [prevls,mnpre]=getvls(avls,pvlindpre);
  
    mintm=17;
    [pstvls,mnpst]=getvls(avls,pvlindpst);
    ax=subplot(427);
    ps.ax=ax
    % ps.axraw=subplot(212)
    % ps.plotraw=1;
    ps(1).ax=ax;
    ps(2).ax=ax;
    ps(1).col='k'
    ps(2).col='r'
    ps(1).ls='-';
    ps(2).ls='-';
    ps(1).wid=2;
    ps(2).wid=2;
%     ps.wnind=[0 0 1 1 1 1 1 1 1 ]
%     ps.wnbnds=[2.37 2.45 2.47 2.47 2.47 2.47 2.47]
%     ps.wnbot=[2.3]
    ps(1).hst=prevls;
    ps(2).hst=pstvls;
    ps(1).edges=avls.edges{1};
    ps(2).edges=avls.edges{1};
    ps(1).plothist=1;
    ps(2).plothist=2;
    ps(1).orient='norm'
    ps(2).orient='norm'
    ps(1).plot_triangles=1;
    ps(2).plot_triangles=2;
    ps(1).mean=mnpre;
    ps(2).mean=mnpst;
    
    plothists2(ps);
end

if(ismember(5,figstoplot))
    pvlindpre=10;
    pvlindpst=13;
    [prevls,mnpre]=getvls(avls,pvlindpre);
  
    mintm=17;
    [pstvls,mnpst]=getvls(avls,pvlindpst);
    ax=subplot(423);
%     ps.ax=ax
    % ps.axraw=subplot(212)
    % ps.plotraw=1;
    ps(1).ax=ax;
    ps(2).ax=ax;
    ps(1).col='k'
    ps(2).col='k'
    ps(1).ls='-';
    ps(2).ls='-';
    ps(1).wid=2;
    ps(2).wid=2;
%     ps.wnind=[0 0 1 1 1 1 1 1 1 ]
%     ps.wnbnds=[2.37 2.45 2.47 2.47 2.47 2.47 2.47]
%     ps.wnbot=[2.3]
    ps(1).hst=prevls;
    ps(2).hst=pstvls;
    ps(1).edges=avls.edges{1};
    ps(2).edges=avls.edges{1};
    ps(1).plothist=1;
    ps(2).plothist=2;
    ps(1).orient='norm'
    ps(2).orient='norm'
    ps(1).plot_triangles=1;
    ps(2).plot_triangles=2;
    ps(1).mean=mnpre;
    ps(2).mean=mnpst;
    
    plothists2(ps);
end





% 
% if(ismember(5,figstoplot))
%     
%     ax=subplot(428);
%     ps.ax=ax
%     % ps.axraw=subplot(212)
%     % ps.plotraw=1;
%     ps.wnax=ax;
%     ps.wnind=[0 0 1 1 1 1 1 1 1 ]
%     ps.wnbnds=[2.37 2.45 2.47 2.47 2.47 2.47 2.47]
%     ps.wnbot=[2.3]
%     plot_tmcourse2revAC(avls,[10 13 16 19 22 24 27 41:44],ps)
% end

function []=plotexspect(exsong,ps)
clear catchvl
clear ax
% 
% figure;
% ax1(1)=subplot(12,1,2:6)
% exsong(1).ax=ax1(1);
% ax1(2)=subplot(12,1,1);

% exsong(2)=exsong(1)
% ax2(1)=subplot(12,1,8:12);
% exsong(2).ax=ax2(1);
% ax2(2)=subplot(12,1,7);


    [dat,fs]=plotcbin2(exsong,ps);
    hold on;
    plot([0 2.5],[6600 6600],'c--')
 
 function [hst,mn]=getvls(avls,pvlind)
    vls=avls.adjvls{1}{pvlind}
    mn=mean(vls(:,2));
    hst=histc(vls(:,2),avls.edges{1});
    hst=hst./sum(hst);
    
