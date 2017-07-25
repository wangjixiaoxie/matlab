%% lt 7/3/17 - given cell array plots each distribution + does anova and/or pairwise comparisons
function lt_plot_MultDist(Yall, xvals, plotanova, plotcol, plotspreadonly, plotanovaonly)


%%
if ~exist('plotanova', 'var')
    plotanova =1;
end

if ~exist('plotcol', 'var')
    plotcol = 'b';
end

if ~exist('plotspreadonly', 'var')
    plotspreadonly=0;
end

if ~exist('xvals', 'var')
    xvals = 1:length(Yall);
end

if ~exist('plotanovaonly', 'var')
    plotanovaonly=0;
end

%%
% --- if any gaps, remove
indstmp = ~cellfun('isempty', Yall);

Yall = Yall(indstmp);
xvals = xvals(indstmp);



%%
% Yall = {Y1, Y2, ...}, where each Y is a vector. will plot in order of
% cell inds

if plotanovaonly==0
if plotspreadonly ==1
%                 sh = plotSpread(ah,data,'xValues',opt.xValues,'xyOri',opt.xyOri);
                sh = plotSpread(Yall, 'xValues', xvals, 'distributionColors', plotcol);
%             set(sh{1},'color',[0,128,255]/255);

else
distributionPlot(Yall, 'xValues', xvals, 'showMM', 4, 'addSpread', 1, 'color', plotcol);
end
end

if plotanova==1
% --- anova *(one way)
Yvec = [];
Group = [];
for i=1:length(Yall)
    if size(Yall{i},1) == 1;
        Yvec = [Yvec; Yall{i}'];
    else
        Yvec = [Yvec; Yall{i}];
    end
    Group = [Group; i*ones(length(Yall{i}),1)];
end

[p] = anovan(Yvec, Group, 'display', 'off');
lt_plot_pvalue(p, 'anova', 2);
end
