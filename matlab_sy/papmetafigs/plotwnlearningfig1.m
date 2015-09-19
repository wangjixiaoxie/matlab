function [combvls]=plotwnlearningfig1();
figure
figstoplot=[4 ]

  avlspath='/oriole6/bk61w42/datasum'
  avlsmat='pathvals1-analdata.mat'
    
  cmd=['cd ' avlspath];
  eval(cmd);
  cmd=['load ' avlsmat];
  eval(cmd);
if(ismember(1,figstoplot))
     
    exsong(1).path=[ avls.pvls{4}]
    exsong(1).cvl=avls.cvl{4}
    exsong(1).fn='bk63w42_011208_0757.3471.cbin';
    exsong(1).bnds=[5 7.1]
    exsong(1).ax=subplot(421)
    plotexspect(exsong);
%     exsong(1).stim_path='/oriole/bk48w74/stimexample'
%     exsong(1).stim_mat='stimexample.mat'
end

if(ismember(3,figstoplot))
     
    exsong(1).path=[ avls.pvls{16}]
    exsong(1).cvl=avls.cvl{16}
    exsong(1).fn='bk63w42_051208_0754.2248.cbin';
    exsong(1).bnds=[3 5.1]
    exsong(1).ax=subplot(425)
    plotexspect(exsong);
%     exsong(1).stim_path='/oriole/bk48w74/stimexample'
%     exsong(1).stim_mat='stimexample.mat'
end


if(ismember(4,figstoplot))
    
    ax=subplot(427);
    ps.ax=ax
    % ps.axraw=subplot(212)
    % ps.plotraw=1;
    axnew=subplot(428);
    ps.wnax=axnew;
    ps.wnind=[0 0 1 1 1 1 1 1 1 ]
    ps.wnbnds=[2.37 2.45 2.47 2.47 2.47 2.47 2.47]
    ps.wnbot=[2.3]
    ps.FILL=1;
    plot_tmcourse2revAC(avls,[1:12],ps)
end

function []=plotexspect(exsong);
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


    [dat,fs]=plotcbin(exsong);
    hold on;
    plot([0 2.5],[6600 6600],'c--')
