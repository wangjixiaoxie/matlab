%% LT 2/1/16 - replot consolidation figure from Warren et al., 2011
lt_figure; hold on;
title('consolidation from Warren et al')

PBSmean=[3.9 3.7 3.5]; % in height in figure
PBSsem=[0.4 0.3 0.3];

MUSCmean=[2.1 2.6 2.92];
MUSCsem=[0.3 0.14 0.2];

    

% normalize all to height of figure at 100%
Norm = 3.68;

PBSmean=PBSmean./Norm;
PBSsem=PBSsem./Norm;

MUSCmean=MUSCmean./Norm;
MUSCsem=MUSCsem./Norm;

X=[1 4 7];
lt_plot_bar(X-0.2, PBSmean, {'Errors',PBSsem, 'BarWidth',0.1})
lt_plot_bar(X+0.2, MUSCmean, {'Errors',MUSCsem, 'Color','k', 'BarWidth',0.1})