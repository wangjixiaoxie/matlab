function lt_seq_dep_pitch_PlotContours(AllDays_RawDatStruct, Params, day);
%% LT 9/7/15 - quick, to plot all syls for a given day, overlaid with time windows used to calculate pitch

count=1;
SubplotsPerFig=9;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

% ====
syls_list=fieldnames(AllDays_RawDatStruct{day}.data);
for i=1:length(syls_list);
    syl=syls_list{i};
    
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title(syl);

PCmat=cell2mat(AllDays_RawDatStruct{day}.data_WithOutlier.(syl)(:,2));

plot(PCmat');

N=size(PCmat,1);

Ylim=ylim;
lt_plot_text(5, Ylim(2)-(Ylim(2)-Ylim(1))/6, ['N = ' num2str(N)])

% put lines for window
twindow=Params.SeqFilter.pc_time_window_list{day}(:,i);

line([twindow(1) twindow(1)], ylim, 'Color','k');
line([twindow(2) twindow(2)], ylim, 'Color','k');

end