function lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT, fieldname_baseneur, fieldname_trainneur)

lt_figure; hold on;
Yall = {};

Yall2 = {}; % {targpre, ssamepre, diffpre, targpost, samepost, diffpost}

minall = [];
maxall = [];

hsplots = [];


plotcols_bird = lt_make_plot_colors(max([DATSTRUCT.AllBirdnum]), 0, 0);

%% pre vs. post (scatter

% -- targ
hsplot = lt_subplot(2,4,1); hold on;
title('targ');
xlabel('base (neural sim)');
ylabel('train');
hsplots = [hsplots hsplot];

inds = DATSTRUCT.AllIsTarg==1;
plotcol = 'k';

X = DATSTRUCT.(fieldname_baseneur)(inds);
Y = DATSTRUCT.(fieldname_trainneur)(inds);
lt_plot_45degScatter(X, Y, plotcol, 1);

Yall{1} = Y-X;
Yall2{1} = X;
Yall2{4} = Y;
% assert(all(~isnan(Y-X) == ~isnan(DATSTRUCT.AllNeurSimChange(inds))));
minall = min([minall X Y]);
maxall = max([maxall X Y]);
% -- use diff color for each bird
birdnums = DATSTRUCT.AllBirdnum(inds);
for i=1:max([DATSTRUCT.AllBirdnum])
   indstmp = birdnums == i;
%    plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i});
   plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i}, 'MarkerFaceColor', plotcols_bird{i});
end


% -- same
hsplot = lt_subplot(2,4,2); hold on;
title('sametype');
xlabel('base');
ylabel('train');
hsplots = [hsplots hsplot];

inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1;
plotcol = 'b';

if sum(inds)>0
X = DATSTRUCT.(fieldname_baseneur)(inds);
Y = DATSTRUCT.(fieldname_trainneur)(inds);
lt_plot_45degScatter(X, Y, plotcol, 1);

Yall{2} = Y-X;
Yall2{2} = X;
Yall2{5} = Y;
minall = min([minall X Y]);
maxall = max([maxall X Y]);
% -- use diff color for each bird
birdnums = DATSTRUCT.AllBirdnum(inds);
for i=1:max([DATSTRUCT.AllBirdnum])
   indstmp = birdnums == i;
%    plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i});
   plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i}, 'MarkerFaceColor', plotcols_bird{i});
end
end


% -- diff
hsplot =lt_subplot(2,4,3); hold on;
title('difftype');
xlabel('base');
ylabel('train');
hsplots = [hsplots hsplot];

inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsDiff==1;
plotcol = 'r';

X = DATSTRUCT.(fieldname_baseneur)(inds);
Y = DATSTRUCT.(fieldname_trainneur)(inds);
lt_plot_45degScatter(X, Y, plotcol, 1);

Yall{3} = Y-X;
Yall2{3} = X;
Yall2{6} = Y;
minall = min([minall X Y]);
maxall = max([maxall X Y]);
% -- use diff color for each bird
birdnums = DATSTRUCT.AllBirdnum(inds);
for i=1:max([DATSTRUCT.AllBirdnum])
   indstmp = birdnums == i;
%    plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i});
   plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i}, 'MarkerFaceColor', plotcols_bird{i});
end

% ---
linkaxes(hsplots, 'xy');
xlim(1.02*[minall maxall]);
ylim(1.02*[minall maxall]);

%% distribtion of pre and post neural similairty

lt_subplot(2,4,4);
title('pre --- post');
ylabel('neural sim')

hold on;
lt_plot_MultDist(Yall2);
% xlim([0 7]);
lt_plot_zeroline;


%% distribution of changes

% ---- all
lt_subplot(2,4,8); hold on;
lt_plot_MultDist(Yall);
lt_plot_zeroline;
ylabel('change in neural sim');

%% change vs. start
hsplots = [];

% -- targ
hsplot = lt_subplot(2,4,5); hold on;
title('targ');
xlabel('base (neural sim)');
ylabel('change in neural sim');
hsplots = [hsplots hsplot];

inds = DATSTRUCT.AllIsTarg==1;
plotcol = 'k';


X = DATSTRUCT.(fieldname_baseneur)(inds);
Y = DATSTRUCT.(fieldname_trainneur)(inds) - DATSTRUCT.(fieldname_baseneur)(inds);
plot(X, Y, 'o', 'Color', plotcol);

lt_plot_zeroline
% -- use diff color for each bird
birdnums = DATSTRUCT.AllBirdnum(inds);
for i=1:max([DATSTRUCT.AllBirdnum])
   indstmp = birdnums == i;
%    plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i});
   plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i}, 'MarkerFaceColor', plotcols_bird{i});
end


% -- same
hsplot = lt_subplot(2,4,6); hold on;
title('same');
xlabel('base (neural sim)');
ylabel('change in neural sim');
hsplots = [hsplots hsplot];

inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1;
plotcol = 'b';

X = DATSTRUCT.(fieldname_baseneur)(inds);
Y = DATSTRUCT.(fieldname_trainneur)(inds) - DATSTRUCT.(fieldname_baseneur)(inds);
plot(X, Y, 'o', 'Color', plotcol);

lt_plot_zeroline
% -- use diff color for each bird
birdnums = DATSTRUCT.AllBirdnum(inds);
for i=1:max([DATSTRUCT.AllBirdnum])
   indstmp = birdnums == i;
%    plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i});
   plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i}, 'MarkerFaceColor', plotcols_bird{i});
end


% -- diff
hsplot = lt_subplot(2,4,7); hold on;
title('diff');
xlabel('base (neural sim)');
ylabel('change in neural sim');
hsplots = [hsplots hsplot];

inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsDiff==1;
plotcol = 'r';

X = DATSTRUCT.(fieldname_baseneur)(inds);
Y = DATSTRUCT.(fieldname_trainneur)(inds) - DATSTRUCT.(fieldname_baseneur)(inds);
plot(X, Y, 'o', 'Color', plotcol);

lt_plot_zeroline
% -- use diff color for each bird
birdnums = DATSTRUCT.AllBirdnum(inds);
for i=1:max([DATSTRUCT.AllBirdnum])
   indstmp = birdnums == i;
%    plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i});
   plot(X(indstmp), Y(indstmp), 'o', 'Color', plotcols_bird{i}, 'MarkerFaceColor', plotcols_bird{i});
end

%--
linkaxes(hsplots, 'xy');
xlim(1.02*[min(DATSTRUCT.(fieldname_baseneur)) max(DATSTRUCT.(fieldname_baseneur))]);
ylim(1.02*[min(DATSTRUCT.(fieldname_trainneur)-DATSTRUCT.(fieldname_baseneur)) ...
    max(DATSTRUCT.(fieldname_trainneur)-DATSTRUCT.(fieldname_baseneur))])


