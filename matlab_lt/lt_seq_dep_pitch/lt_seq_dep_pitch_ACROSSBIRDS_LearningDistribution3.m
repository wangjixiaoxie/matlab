function lt_seq_dep_pitch_ACROSSBIRDS_LearningDistribution3(SeqDepPitch_AcrossBirds, PARAMS)


%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);



%% Initiate plots
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


%% =========== plot line for each expt (targ -- same -- diff).
% overlay with  bars + baseline drift distribution

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            Learn_targ_all=[];
            Learn_same_all=[];
            Learn_diff_all=[];


for i=1:NumBirds
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        % color for this expt
        plotcol=[rand rand rand];
        jitter=rand*0.2;
        
        % --------------- SANITY CHECK --------------------------
        % --- confirm that the previously extracted syl fields are correct
        list1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTDS;
        list2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTSS;
        list3=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS;
        list4=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS;
        sylsunique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        newlist=[list1 list2 list3 list4];
        newlist=unique(newlist);
        newlist=[newlist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl]; % this should be identical to syls unqiue
        
        assert(length(newlist)==length(sylsunique), 'PROBLEM');
        assert(length(intersect(newlist, sylsunique))==length(sylsunique), 'PROBLEM'); % conf irms that intersect is same length as either cell array
        % -----------------------------------------
        
        % --- learn dir
        targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % ====================== plot targ, sametype, diff type
        % - targ
        x=2;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        learn_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        learn_targ=learn_targ*targdir;
        
%         plot(x+0.2, learn_targ, 'ok');
        plot(x+jitter, learn_targ, 'ok');
        
        Learn_targ_all=[Learn_targ_all learn_targ];
        
        % - same
        x=1;
        
        sylsSame=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
        
        for j=1:length(sylsSame)
            syl=sylsSame{j};
            
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learn=learn*targdir;
            
            % -- plot line
            %             plot([x+0.2 2.2], [learn learn_targ], '-', 'Color', plotcol);
            plot([x+jitter 2+jitter], [learn learn_targ], '-', 'Color', [0.8 0.8 0.8]);
            
            % -- plot dot
%             plot(x+0.2, learn, 'ob');
            plot(x+jitter, learn, 'ob');
            
            % --- collect values to get mean
            Learn_same_all=[Learn_same_all learn];
            
        end
        
        
        % - diff
        x=3;
        
        sylsDiff=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTSS];
        
        for j=1:length(sylsDiff)
            syl=sylsDiff{j};
            
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learn=learn*targdir;
            
            % -- plot line
            %             plot([x+0.2 2.2], [learn learn_targ], '-', 'Color', plotcol);
%             plot([x+0.2 2.2], [learn learn_targ], '-', 'Color', [0.8 0.8 0.8]);
            plot([x+jitter 2+jitter], [learn learn_targ], '-', 'Color', [0.8 0.8 0.8]);
            
            % -- plot dot
%             plot(x+0.2, learn, 'or');
            plot(x+jitter, learn, 'or');

            Learn_diff_all=[Learn_diff_all learn];
        end
        
                
    end
end

% ==== plot mean across all syls
lt_plot_bar(1, mean(Learn_same_all), {'Errors', lt_sem(Learn_same_all), 'Color','b'}); 
lt_plot_bar(2, mean(Learn_targ_all), {'Errors', lt_sem(Learn_targ_all), 'Color','k'}); 
lt_plot_bar(3, mean(Learn_diff_all), {'Errors', lt_sem(Learn_diff_all), 'Color','r'}); 


% === plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');


%% =========== plot line for each expt (targ -- same -- diff). [TRY DIFFERENT ORDER IN PLOT]
% overlay with  bars + baseline drift distribution

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            Learn_targ_all=[];
            Learn_same_all=[];
            Learn_diff_all=[];


for i=1:NumBirds
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        % color for this expt
        plotcol=[rand rand rand];
        jitter=rand*0.2;
        
        % --------------- SANITY CHECK --------------------------
        % --- confirm that the previously extracted syl fields are correct
        list1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTDS;
        list2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTSS;
        list3=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS;
        list4=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS;
        sylsunique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        newlist=[list1 list2 list3 list4];
        newlist=unique(newlist);
        newlist=[newlist SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl]; % this should be identical to syls unqiue
        
        assert(length(newlist)==length(sylsunique), 'PROBLEM');
        assert(length(intersect(newlist, sylsunique))==length(sylsunique), 'PROBLEM'); % conf irms that intersect is same length as either cell array
        % -----------------------------------------
        
        % --- learn dir
        targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % ====================== plot targ, sametype, diff type
        % - targ
        x=1;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        learn_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        learn_targ=learn_targ*targdir;
        
%         plot(x+0.2, learn_targ, 'ok');
        plot(x+jitter, learn_targ, 'o', 'Color', 'k');
        
        Learn_targ_all=[Learn_targ_all learn_targ];
        
        % - same
        x=2;
        
        sylsSame=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
        
        for j=1:length(sylsSame)
            syl=sylsSame{j};
            
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learn=learn*targdir;
            
            % -- plot line
            %             plot([x+0.2 2.2], [learn learn_targ], '-', 'Color', plotcol);
%             plot([x+jitter 1+jitter], [learn learn_targ], '-', 'Color', [0.8 0.8 0.8]);
%             plot([x+jitter 1+jitter], [learn learn_targ], '-', 'Color', plotcol);
            
            % -- plot dot
%             plot(x+0.2, learn, 'ob');
            plot(x+jitter, learn, 'o', 'Color','b');
            
            % --- collect values to get mean
            Learn_same_all=[Learn_same_all learn];
            
        end
        
        
        % - diff
        x=3;
        
        sylsDiff=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTSS];
        
        for j=1:length(sylsDiff)
            syl=sylsDiff{j};
            
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learn=learn*targdir;
            
            % -- plot line
            %             plot([x+0.2 2.2], [learn learn_targ], '-', 'Color', plotcol);
%             plot([x+0.2 1.2], [learn learn_targ], '-', 'Color', [0.8 0.8 0.8]);
%             plot([x+jitter 1+jitter], [learn learn_targ], '-', 'Color', plotcol);
            
            % -- plot dot
%             plot(x+0.2, learn, 'or');
            plot(x+jitter, learn, 'o', 'Color', 'r');

            Learn_diff_all=[Learn_diff_all learn];
        end
        
                
    end
end

% ==== plot mean across all syls
lt_plot_bar(2, mean(Learn_same_all), {'Errors', lt_sem(Learn_same_all), 'Color','b'}); 
lt_plot_bar(1, mean(Learn_targ_all), {'Errors', lt_sem(Learn_targ_all), 'Color','k'}); 
lt_plot_bar(3, mean(Learn_diff_all), {'Errors', lt_sem(Learn_diff_all), 'Color','r'}); 


% === plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');


%% ================= DISTRIBUTION PLOTS, BOX AND WHISKER, ETC

            Learn_targ_all=[];
            Learn_same_all=[];
            Learn_diff_all=[];


for i=1:NumBirds
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        % color for this expt
%         plotcol=[rand rand rand];
        jitter=rand*0.2;

        % --- learn dir
        targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % ====================== plot targ, sametype, diff type
        % - targ
        x=1;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        learn_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        learn_targ=learn_targ*targdir;
                
        Learn_targ_all=[Learn_targ_all learn_targ];
        
        % - same
        x=2;
        
        sylsSame=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
        
        for j=1:length(sylsSame)
            syl=sylsSame{j};
            
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learn=learn*targdir;
            
            % --- collect values to get mean
            Learn_same_all=[Learn_same_all learn];
        end
        
        
        % - diff
        x=3;
        
        sylsDiff=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTSS];
        
        for j=1:length(sylsDiff)
            syl=sylsDiff{j};
            
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learn=learn*targdir;
            
            Learn_diff_all=[Learn_diff_all learn];
        end
    end
end


% ======================== BOX AND WHISKER
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% ==== plot mean across all syls
X=[Learn_targ_all, Learn_same_all, Learn_diff_all];
G=[ones(1,length(Learn_targ_all)), 2*ones(1,length(Learn_same_all)), 3*ones(1,length(Learn_diff_all))];
boxplot(X, G, 'color','kbr', 'symbol', '.', 'jitter', 0.1); 


%  plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');



% ============================= DISTRIBUTION PLOT (ksd)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

distributionPlot({Learn_targ_all', Learn_same_all', Learn_diff_all'}, 'color', {'k', 'b', 'r'})

%  plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');


% ============================= DISTRIBUTION PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

distributionPlot({Learn_targ_all', Learn_same_all', Learn_diff_all'}, 'color', {'k', 'b', 'r'}, 'histOpt', 0)

%  plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');


% ============================ PLOT SCATTER WITH NOT HIST
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

handles=plotSpread({Learn_targ_all', Learn_same_all', Learn_diff_all'}, 'distributionColors', {'k', 'b', 'r'}, 'showMM', 4);

set(handles{2}(1), 'Color','k');
set(handles{2}(2), 'Color','k');


%  plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');

line(xlim, [0 0], 'Color' ,'k');

% -- sample sizes;
n1=length(Learn_targ_all);
n2=length(Learn_same_all);
n3=length(Learn_diff_all);

lt_plot_annotation(1, ['ntarg=' num2str(n1) '; nsame=' num2str(n2) '; ndiff=' num2str(n3)]);


% --- stats, directly compare using rank sum
p1v2=ranksum(Learn_targ_all, Learn_same_all)
p2v3=ranksum(Learn_same_all, Learn_diff_all)
p1v3=ranksum(Learn_targ_all, Learn_diff_all)

p1=signrank(Learn_targ_all)
p2=signrank(Learn_same_all)
p3=signrank(Learn_diff_all)

% test variance vs. baseline drift
baselineDrift=PARAMS.baselineDrift.AllDriftVals;
[h, p] = kstest2(baselineDrift, Learn_diff_all)





% ============================ OVERLAY DISTRIBUTIONS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

Xcenters=linspace(min([Learn_targ_all, Learn_same_all, Learn_diff_all]), max([Learn_targ_all, Learn_same_all, Learn_diff_all]),33);

lt_plot_histogram(Learn_targ_all, Xcenters, 1, 1, '', 1, 'k');
lt_plot_histogram(Learn_same_all, Xcenters, 1, 1, '', 1, 'b');
lt_plot_histogram(Learn_diff_all, Xcenters, 1, 1, '', 1, 'r');
% 
% plotSpread({Learn_targ_all', Learn_same_all', Learn_diff_all'}, 'distributionColors', {'k', 'b', 'r'}, 'showMM', 4)

%  plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');


%% scatter of distribution as before, but showing gaussian of baseline drift
lt_figure; hold on;

% ============================ PLOT SCATTER WITH NOT HIST

handles=plotSpread({Learn_targ_all', Learn_same_all', Learn_diff_all'}, 'distributionColors', {'k', 'b', 'r'}, 'showMM', 4);

set(handles{2}(1), 'Color','k');
set(handles{2}(2), 'Color','k');


% === fit gaussian to baseline drift vals
pd=fitdist(PARAMS.baselineDrift.AllDriftVals', 'Normal');

% get cdf
xtmp=-100:0.1:100;
Ycdf=cdf(pd, xtmp);
lt_subplot(2,2,1); subtitle('cdf of baseline drift, fit to gaussian');
plot(xtmp, Ycdf);

% -- get val for 2.5 and 97.5 percentiles
prctile1=2.5;
[C, ind]=min(abs(Ycdf-prctile1/100));
lowerBound=xtmp(ind);

prctile2=97.5;
[C, ind]=min(abs(Ycdf-prctile2/100));
upperBound=xtmp(ind);



%  plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');

line(xlim, [0 0], 'Color' ,'k');

% -- sample sizes;
n1=length(Learn_targ_all);
n2=length(Learn_same_all);
n3=length(Learn_diff_all);

lt_plot_annotation(1, ['ntarg=' num2str(n1) '; nsame=' num2str(n2) '; ndiff=' num2str(n3)]);


%% PLOT ONE POINT FOR EACH EXPERIMENT, THEN PUT LINE
% --- ONLY PLOT IF EXPT HAS BOTH SAME AND DIFF TYPE

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
for i=1:NumBirds
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        % color for this expt
        plotcol=[rand rand rand];
        jitter=rand*0.2;
        
        
        % --- learn dir
        targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        % ====================== plot targ, sametype, diff type
        % - targ
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        learn_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        learn_targ=learn_targ*targdir;
        
        
        % - same
        sylsSame=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
        
        Learn_same_all=[];
        
        for j=1:length(sylsSame)
            syl=sylsSame{j};
            
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learn=learn*targdir;

            % --- collect values to get mean
            Learn_same_all=[Learn_same_all learn];
        end
        learn_same_mean=mean(Learn_same_all);
        
        
        
        % - diff
        sylsDiff=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTSS];
        
        Learn_diff_all=[];
        
        for j=1:length(sylsDiff)
            syl=sylsDiff{j};
            
            learn=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learn=learn*targdir;
            
            Learn_diff_all=[Learn_diff_all learn];
        end
        
        learn_diff_mean=mean(Learn_diff_all);
        
        % ====== PLOT FOR THIS EXPT
        % - only plot if expt has both same and diff type
        if isnan(learn_same_mean) | isnan(learn_diff_mean)
        else
            
            plot([1 2 3], [learn_targ learn_same_mean learn_diff_mean], 'o-');
            
        end
        
                
    end
end

% % ==== plot mean across all syls
% lt_plot_bar(1, mean(Learn_same_all), {'Errors', lt_sem(Learn_same_all), 'Color','b'}); 
% lt_plot_bar(2, mean(Learn_targ_all), {'Errors', lt_sem(Learn_targ_all), 'Color','k'}); 
% lt_plot_bar(3, mean(Learn_diff_all), {'Errors', lt_sem(Learn_diff_all), 'Color','r'}); 


% === plot lines for baseline drift
lowerBound=PARAMS.baselineDrift.prctiles_2_5and97_5(1);
upperBound=PARAMS.baselineDrift.prctiles_2_5and97_5(2);

line(xlim, [lowerBound lowerBound], 'Color','k','LineStyle','--');
line(xlim, [upperBound upperBound], 'Color','k','LineStyle','--');

title('only expt with both same and diff, average over those classes');

