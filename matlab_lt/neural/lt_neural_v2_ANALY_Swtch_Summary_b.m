function lt_neural_v2_ANALY_Swtch_Summary_a(DATSTRUCT, fieldname_baseneur, fieldname_trainneur)

lt_figure; hold on;
Yall = {};

Yall2 = {}; % {targpre, ssamepre, diffpre, targpost, samepost, diffpost}

minall = [];
maxall = [];

hsplots = [];

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
assert(all(~isnan(Y-X) == ~isnan(DATSTRUCT.AllNeurSimChange(inds))));
minall = min([minall X Y]);
maxall = max([maxall X Y]);


% -- nontarg
hsplot = lt_subplot(2,4,2); hold on;
title('nontarg');
xlabel('base');
ylabel('train');
hsplots = [hsplots hsplot];

inds = DATSTRUCT.AllIsTarg==0;
plotcol = 'm';

X = DATSTRUCT.(fieldname_baseneur)(inds);
Y = DATSTRUCT.(fieldname_trainneur)(inds);
lt_plot_45degScatter(X, Y, plotcol, 1);

Yall{2} = Y-X;
Yall2{2} = X;
Yall2{5} = Y;
minall = min([minall X Y]);
maxall = max([maxall X Y]);

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


% -- nontarg
hsplot = lt_subplot(2,4,6); hold on;
title('nontarg');
xlabel('base (neural sim)');
ylabel('change in neural sim');
hsplots = [hsplots hsplot];

inds = DATSTRUCT.AllIsTarg==0;
plotcol = 'm';

X = DATSTRUCT.(fieldname_baseneur)(inds);
Y = DATSTRUCT.(fieldname_trainneur)(inds) - DATSTRUCT.(fieldname_baseneur)(inds);
plot(X, Y, 'o', 'Color', plotcol);

lt_plot_zeroline


%--
linkaxes(hsplots, 'xy');
xlim(1.02*[min(DATSTRUCT.(fieldname_baseneur)) max(DATSTRUCT.(fieldname_baseneur))]);
ylim(1.02*[min(DATSTRUCT.(fieldname_trainneur)-DATSTRUCT.(fieldname_baseneur)) ...
    max(DATSTRUCT.(fieldname_trainneur)-DATSTRUCT.(fieldname_baseneur))])

