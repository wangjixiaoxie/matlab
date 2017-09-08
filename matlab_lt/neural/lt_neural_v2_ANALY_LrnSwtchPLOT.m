function lt_neural_v2_ANALY_LrnSwtchPLOT(MOTIFSTATS_Compiled, SwitchStruct)

% NOTE: currently taking 20 renditions at the end of training 
Nrends = 20;
minrends = 3; % pre and post
plotRaw =1; % for each swi5tch

%%
assert(length(MOTIFSTATS_Compiled.birds) == length(SwitchStruct.bird), 'numbirds dont match')


%% 


figcount=1;
subplotrows=3;
subplotcols=5;
fignums_alreadyused=[];
hfigs=[];


numbirds = length(SwitchStruct.bird);

LearnAll = []; % switches x [targ, same, diff];

for i=1:numbirds
    birdname = SwitchStruct.bird(i).birdname;
    numexpts = length(SwitchStruct.bird(i).exptnum);
    
    for ii=1:numexpts
       exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
       
        for iii=1:length(SwitchStruct.bird(i).exptnum(ii).switchlist)
           
            if isempty(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs)
                keyboard
            end
            
            % --- skip if targs diff
            if SwitchStruct.bird(i).exptnum(ii).switchlist(iii).targsAreSameSyl==0
                continue
            end
            
            % --- skip if bidir learning
            if length(unique(cell2mat(SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs(2:2:end))))>1
                % then bidir (or one targ no change)
                continue
            end
           
            % --- skip if no data (neurons)
            if ~any([SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspostsongs] & ...
                    [SwitchStruct.bird(i).exptnum(ii).switchlist(iii).neuron.haspresongs])
                continue
            end
            
            % --- skip if no switch
            learndir = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).learningDirs{2};
            if learndir==0
            continue
            end
            
            
            % ------------
            swdnum_pre = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_previous;
            swdnum_post = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum_next;
            swdnum_this = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).switchdnum;
            
            
            % ===================== COLLECT FFVALS (COLLAPSE WITHIN
            % CATEGORIES BY TAKING MEDIAN);
            motifstats = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS;
            motiflist = motifstats.params.motif_regexpr_str;
            
            Learning_targ = [];
            Learning_sametype = [];
            Learning_difftype = [];
            for j=1:length(motiflist)
                syl = motiflist{j};
                
                indtmp = strcmp({SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl.sylname}, ...
                    syl);
                
                assert(j == find(indtmp), 'sasadf');
                
                tvals = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(indtmp).tvals;
                ffvals = SwitchStruct.bird(i).exptnum(ii).switchlist(iii).STATS_motifsyl(indtmp).ffvals;
                
                % --- get diff from base
                baseInds = find(tvals>swdnum_pre & tvals<swdnum_this);
                trainInds = find(tvals>swdnum_this & tvals<swdnum_post);
                
%                 if any(strcmp(motifstats.params.TargSyls, syl))
%                     keyboard 
%                 end
%                 
                if length(baseInds)<minrends | length(trainInds)<minrends
                    continue
                end
                
                ffbasemean = mean(ffvals(baseInds));
                ffbaseSD = std(ffvals(baseInds));
                
                if isnan(ffbasemean)
                    continue
                end
                
                if length(trainInds)>=Nrends
                ffvals_train = ffvals(trainInds(end-Nrends+1:end));
                else
                    ffvals_train = ffvals(trainInds(1:end));
                end
                 ffvals_train = (ffvals_train - ffbasemean)./ffbaseSD;
               
                % ---- final
                learndiff = mean(ffvals_train);
                learndiff = learndiff*learndir;
                
                % --- output
                if any(strcmp(motifstats.params.TargSyls, syl))
                    Learning_targ = [Learning_targ learndiff];
                elseif any(strcmp(motifstats.params.SameTypeSyls, syl))
                    Learning_sametype = [Learning_sametype learndiff];
                elseif any(strcmp(motifstats.params.DiffTypeSyls, syl))
                    Learning_difftype = [Learning_difftype learndiff];
                end                    
            end
            
            % ==================== PLOT
            if plotRaw ==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-sw' num2str(iii)])
                
                
                plot(1, Learning_targ, 'ok');
                if ~isempty(Learning_sametype)
                plot(2, Learning_sametype, 'ob');
                end
                if ~isempty(Learning_difftype)
                plot(3, Learning_difftype, 'or');
                end
                xlim([0 4]);
                lt_plot_zeroline;
                         
            end
           
            % =============================== GET MEDIAN FOR EACH CATEGORY
%             if any(isnan(Learning_difftype))
%                 disp(Learning_difftype)
%             end
%             LearnAll = [LearnAll; [median(Learning_targ) median(Learning_sametype) median(Learning_difftype)]];
            LearnAll = [LearnAll; [mean(Learning_targ) mean(Learning_sametype) mean(Learning_difftype)]];
            
            
            
        end
    end
end

%% ==== PLOT - one line for each switch (medians)

lt_figure; hold on;
title('each line one switch (mean within class)');
ylabel('learning (z)');

plot(1:3, LearnAll', '-b');
% lt_plot(1.2:3.2, nanmean(LearnAll, 1), {'Errors', lt_sem(LearnAll), 'LineStyle', '-', 'Color', 'k'});
lt_plot_bar(1.2:3.2, nanmean(LearnAll, 1), {'Errors', lt_sem(LearnAll), 'LineStyle', '-', 'Color', 'k'});
lt_plot_zeroline;
xlim([0 4]);

% -- test significance
for j=1:size(LearnAll,2)
   p = signrank(LearnAll(:,j)); 
   lt_plot_text(j, 1.1*max(LearnAll(:,j)), ['p=' num2str(p)], 'r')
    
end

set(gca, 'XTick', 1:3, 'XTickLabel', {'targ', 'same', 'diff'});





