%% how to run:

% --- OTUPUTS
% 

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Histograms: pitch diff (raw) (all nontargs)');
grid on;


% get histogram centers
[~, Xcenters]=lt_plot_histogram([Y(:,1)' Y(:,2)'], '',0,0);
xNames={'baseline','learning'};
distributionPlot(Y,'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);


% ======= 9b) HISTOGRAM (in direction of learning for that experiment)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
xlabel('pitch diff (flipped if targ learn is neg)');
grid on;

targ_learn_sign=TargLearnDir_All(inds);

Y_flipped=Y.*repmat(targ_learn_sign',1,2);

% get histogram centers
xNames={'baseline','learning'};
distributionPlot(Y_flipped,'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);



% ======= 9b) 
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
xlabel('pitch (in targ dir; base sign matched to learning)');
grid on;

% -- match sign of baseline to learning
Y_learn_signs=sign(Y_flipped(:,2));
Y_base=Y_flipped(:,1);
Y_base_signmatch=abs(Y_base).*Y_learn_signs;

Y_flipped_BaseSignMatched=[Y_base_signmatch Y_flipped(:,2)];

% get histogram centers
xNames={'baseline','learning'};
distributionPlot(Y_flipped_BaseSignMatched,'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);





% ======= 9c) POSITIVE GEN
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('POSITIVE GENERALIZERS');
xlabel('pitch diff (in targ dir; base sign matched to learning)');
grid on;

inds_tmp=Y_flipped_BaseSignMatched(:,2)>0;


% get histogram centers
xNames={'baseline','learning'};
distributionPlot(Y_flipped_BaseSignMatched(inds_tmp,:),'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);



% ======= 9c) NEGATIVE GEN
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('NEGATIVE GENERALIZERS');
xlabel('pitch diff (in targ dir; base sign matched to learning)');
grid on;

inds_tmp=Y_flipped_BaseSignMatched(:,2)<0;


% get histogram centers
if any(inds_tmp)
xNames={'baseline','learning'};
try
distributionPlot(Y_flipped_BaseSignMatched(inds_tmp,:),'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor', Xcenters);
catch err
 distributionPlot(Y_flipped_BaseSignMatched(inds_tmp,:),'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames);
end
end

% ======= 9c) HISTOGRAM (abs)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
xlabel('pitch diff (abs)');
grid on;

Y_abs=abs(Y);

% get histogram centers
xNames={'baseline','learning'};
distributionPlot(Y_abs,'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);





%% CDFs
% ======= 9c) CDF
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs] (dash: bline, solid: learn)');
xlabel('pitch diff (raw)');
grid on;

% -- RAW
% baseline
h=lt_plot_cdf(Y(:,1),'k');
set(h,'LineStyle','--')
% learning (raw)
lt_plot_cdf(Y(:,2),'k');

% --- KS test
[~, p] = kstest2(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'kstest');


% ======= 9c) CDF
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
xlabel('pitch diff (flipped if targ learn neg)');
grid on;

% -- FLIPPED
% baseline
h=lt_plot_cdf(Y_flipped(:,1),'k');
set(h,'LineStyle','--')
% learning
lt_plot_cdf(Y_flipped(:,2),'k');

% --- KS test
[~, p] = kstest2(Y_flipped(:,1), Y_flipped(:,2));
lt_plot_pvalue(p, 'kstest');


% ======= 9c) CDF
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
xlabel('pitch (in targ dir; base sign matched to learning)');
grid on;

% -- FLIPPED
% baseline
h=lt_plot_cdf(Y_flipped_BaseSignMatched(:,1),'k');
set(h,'LineStyle','--')
% learning
lt_plot_cdf(Y_flipped_BaseSignMatched(:,2),'k');

% --- KS test
[~, p] = kstest2(Y_flipped_BaseSignMatched(:,1), Y_flipped_BaseSignMatched(:,2));
lt_plot_pvalue(p, 'kstest');



% ======= 9c) CDF
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
xlabel('pitch diff (abs)');
grid on;

% -- ABS
% baseline
h=lt_plot_cdf(Y_abs(:,1),'k');
set(h,'LineStyle','--')
% learning (raw)
lt_plot_cdf(Y_abs(:,2),'k');

% --- KS test
[~, p] = kstest2(Y_abs(:,1), Y_abs(:,2));
lt_plot_pvalue(p, 'kstest');


%% PAIRED DIFFS
% ======= 9c) CDF (paired diff)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
grid on;

% --- raw
xlabel('diff raw pitch diff (learning - base)');
Y_diffpaired=Y(:,2)-Y(:,1);
lt_plot_cdf(Y_diffpaired,'k');
p=signrank(Y(:,2), Y(:,1));
lt_plot_pvalue(p, 'two-sided sign rank(raw)');




% ======= 9c) CDF (paired diff)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
grid on;

% --- flipped
xlabel('diff pitch diff (flipped if learn neg) (learning - base)');
Y_diffpaired_flipped=Y_flipped(:,2)-Y_flipped(:,1);
lt_plot_cdf(Y_diffpaired_flipped,'r');

% sign rank
p=signrank(Y_flipped(:,2), Y_flipped(:,1), 'tail','right');
lt_plot_pvalue(p, 'one-sided (learn>base) sign rank', 1);

% sign rank
p=signrank(Y_flipped(:,2), Y_flipped(:,1), 'tail','left');
lt_plot_pvalue(p, 'one-sided (learn<base) sign rank', 2);


% --- OUTPUT EFFECT SIZE
OUTPUT.Paired_analyses.flipped_if_targ_neg.ff_diff_mean=mean(Y_flipped(:,2)-Y_flipped(:,1));
OUTPUT.Paired_analyses.flipped_if_targ_neg.ff_diff_sem=lt_sem(Y_flipped(:,2)-Y_flipped(:,1));
p=signrank(Y_flipped(:,2), Y_flipped(:,1), 'tail','right');
OUTPUT.Paired_analyses.flipped_if_targ_neg.ff_diff_p_learn_greater=p;
p=signrank(Y_flipped(:,2), Y_flipped(:,1), 'tail','left');
OUTPUT.Paired_analyses.flipped_if_targ_neg.ff_diff_p_learn_smaller=p;
p=signrank(Y_flipped(:,2), Y_flipped(:,1));
OUTPUT.Paired_analyses.flipped_if_targ_neg.ff_diff_p_notail=p;




% ======= 9c) CDF (paired diff)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
xlabel('pitch (in targ dir; base sign matched to learning)');
grid on;

% --- flipped
Y_diffpaired_flipped_BaseSignMatched=Y_flipped_BaseSignMatched(:,2)-Y_flipped_BaseSignMatched(:,1);
lt_plot_cdf(Y_diffpaired_flipped_BaseSignMatched,'c');


% sign rank
p=signrank(Y_flipped_BaseSignMatched(:,2), Y_flipped_BaseSignMatched(:,1), 'tail','right');
lt_plot_pvalue(p, 'one-sided (learn>base) sign rank', 1);

% sign rank
p=signrank(Y_flipped_BaseSignMatched(:,2), Y_flipped_BaseSignMatched(:,1), 'tail','left');
lt_plot_pvalue(p, 'one-sided (learn<base) sign rank', 2);




% ======= 9c) CDF (paired diff)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('POSITIVE GENERALIZERS');
xlabel('pitch diff (in targ dir; base sign matched to learning)');
grid on;

% --- flipped
Y_diffpaired_flipped_BaseSignMatched=Y_flipped_BaseSignMatched(:,2)-Y_flipped_BaseSignMatched(:,1);

inds_tmp=Y_flipped_BaseSignMatched(:,2)>0;
Z=Y_diffpaired_flipped_BaseSignMatched(inds_tmp,:);

lt_plot_cdf(Z,'c');

% sign rank
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1), 'tail','right');
lt_plot_pvalue(p, 'one-sided (learn>base) sign rank', 1);

% sign rank
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1), 'tail','left');
lt_plot_pvalue(p, 'one-sided (learn<base) sign rank', 2);

% --- OUTPUT EFFECT SIZE
OUTPUT.Paired_analyses.positive_vs_basematched.ff_diff_mean=mean(Z);
OUTPUT.Paired_analyses.positive_vs_basematched.ff_diff_sem=lt_sem(Z);
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1), 'tail','right');
OUTPUT.Paired_analyses.positive_vs_basematched.ff_diff_p_learn_greater=p;
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1), 'tail','left');
OUTPUT.Paired_analyses.positive_vs_basematched.ff_diff_p_learn_smaller=p;
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1));
OUTPUT.Paired_analyses.positive_vs_basematched.ff_diff_p_notail=p;



% ======= 9c) CDF (paired diff)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('NEGATIVE GENERALIZERS');
xlabel('pitch diff (in targ dir; base sign matched to learning)');
grid on;

% --- flipped
Y_diffpaired_flipped_BaseSignMatched=Y_flipped_BaseSignMatched(:,2)-Y_flipped_BaseSignMatched(:,1);

inds_tmp=Y_flipped_BaseSignMatched(:,2)<0;
Z=Y_diffpaired_flipped_BaseSignMatched(inds_tmp,:);

if any(inds_tmp)
lt_plot_cdf(Z,'c');

% sign rank
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1), 'tail','right');
lt_plot_pvalue(p, 'one-sided (learn>base) sign rank', 1);

% sign rank
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1), 'tail','left');
lt_plot_pvalue(p, 'one-sided (learn<base) sign rank', 2);

% --- OUTPUT EFFECT SIZE
OUTPUT.Paired_analyses.negative_vs_basematched.ff_diff_mean=mean(Z);
OUTPUT.Paired_analyses.negative_vs_basematched.ff_diff_sem=lt_sem(Z);
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1), 'tail','right');
OUTPUT.Paired_analyses.negative_vs_basematched.ff_diff_p_learn_greater=p;
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1), 'tail','left');
OUTPUT.Paired_analyses.negative_vs_basematched.ff_diff_p_learn_smaller=p;
p=signrank(Y_flipped_BaseSignMatched(inds_tmp,2), Y_flipped_BaseSignMatched(inds_tmp,1));
OUTPUT.Paired_analyses.negative_vs_basematched.ff_diff_p_notail=p;
end

% ======= 9c) CDF (paired diff)
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('[nontargs]');
grid on;

% --- abs
xlabel('diff pitch diff (abs) (learning - base)');
Y_diffpaired_abs=Y_abs(:,2)-Y_abs(:,1);
h=lt_plot_cdf(Y_diffpaired_abs,'b');

% sign rank
p=signrank(Y_abs(:,2), Y_abs(:,1), 'tail','right');
lt_plot_pvalue(p, 'one-sided (learn>base) sign rank', 1);

% sign rank
p=signrank(Y_abs(:,2), Y_abs(:,1), 'tail','left');
lt_plot_pvalue(p, 'one-sided (learn<base) sign rank', 2);

% --- OUTPUT EFFECT SIZE
OUTPUT.Paired_analyses.abs_val.ff_diff_mean=mean(Y_abs(:,2)-Y_abs(:,1));
OUTPUT.Paired_analyses.abs_val.ff_diff_sem=lt_sem(Y_abs(:,2)-Y_abs(:,1));
p=signrank(Y_abs(:,2), Y_abs(:,1), 'tail','right');
OUTPUT.Paired_analyses.abs_val.ff_diff_p_learn_greater=p;
p=signrank(Y_abs(:,2), Y_abs(:,1), 'tail','left');
OUTPUT.Paired_analyses.abs_val.ff_diff_p_learn_smaller=p;
p=signrank(Y_abs(:,2), Y_abs(:,1));
OUTPUT.Paired_analyses.abs_val.ff_diff_p_notail=p;


%% ======= 9) HISTOGRAM [OLD]
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('Histograms: pitch diff (raw) (similar)');
% grid on;
% 
% % similar
% inds=Target_All==0 & Similar_All==1;
% 
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% 
% % get histogram centers
% [~, Xcenters]=lt_plot_histogram([FFdiff_base_All(inds) FFdiff_learn_All(inds)], '',0,0);
% xNames={'baseline','learning'};
% distributionPlot(Y,'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);
% 
% 
% % ======= 9) HISTOGRAM
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('Histograms: pitch diff (raw) (different)');
% grid on;
% 
% % similar
% inds=Target_All==0 & Similar_All==0;
% 
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% 
% % get histogram centers
% [~, Xcenters]=lt_plot_histogram([FFdiff_base_All(inds) FFdiff_learn_All(inds)], '',0,0);
% xNames={'baseline','learning'};
% distributionPlot(Y,'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);
% 
% 
% % ======= 9) HISTOGRAM
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('Histograms: pitch diff (raw) (presyl similar)');
% grid on;
% 
% % inds
% inds=Target_All==0 & PreSimilar_All==1;
% 
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% 
% % get histogram centers
% [~, Xcenters]=lt_plot_histogram([FFdiff_base_All(inds) FFdiff_learn_All(inds)], '',0,0);
% xNames={'baseline','learning'};
% distributionPlot(Y,'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);
% 
% % ======= 9) HISTOGRAM
% [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
% title('Histograms: pitch diff (raw) (presyl diff)');
% grid on;
% 
% % inds
% inds=Target_All==0 & PreSimilar_All==0;
% 
% Y=[FFdiff_base_All(inds)' FFdiff_learn_All(inds)'];
% 
% % get histogram centers
% [~, Xcenters]=lt_plot_histogram([FFdiff_base_All(inds) FFdiff_learn_All(inds)], '',0,0);
% xNames={'baseline','learning'};
% distributionPlot(Y,'histOpt',0, 'showMM', 6, 'xyOri','flipped', 'histOri','right','xNames',xNames,'divFactor',Xcenters);
% 
% 
% 
% 