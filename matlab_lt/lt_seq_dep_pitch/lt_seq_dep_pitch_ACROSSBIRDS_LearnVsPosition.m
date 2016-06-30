function lt_seq_dep_pitch_ACROSSBIRDS_LearnVsPosition(SeqDepPitch_AcrossBirds, Params, norm_all_to_targ, plot_error_bars_on_raw, equal_y_axis)


% === REMOVED, becuase I do this at the start of the seq dep pitch code.
% %% SORT OUT ONLY THE EXPERIMENTS WITH ONLY ONE TARG, OR STARTED OUT WITH ONLY ONE TARG
% 
% 
% NumBirds=length(SeqDepPitch_AcrossBirds.birds);
% 
% % remove experiments - first figure out what inds to remove
% expts_to_remove=[];
% for i=1:NumBirds;
%     NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%     
%     expts_to_remove{i}=[];
%     for ii=1:NumExperiments;
%         numtargs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs;
%         
%         if numtargs>1;
%             % mult targets. note down index and remove later
%             expts_to_remove{i}=[expts_to_remove{i} ii];
%         end
%     end
% end
% 
% % actually remove them
% for i=1:length(expts_to_remove);
%     if ~isempty(expts_to_remove{i});
%         SeqDepPitch_AcrossBirds.birds{i}.experiment(expts_to_remove{i})=[];
%         disp(['TO GET ONLY ONE TARGETS: removed: bird: ' num2str(i) '; expt: ' num2str(expts_to_remove{i})]);
%     end
% end
% 
% 
% PARAMS.global.did_take_only_onetarg_expts=1;
%% ===
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% PLOT LEARNING REL TO DISTANCE IN MOTIF
write_syl=0;

lt_figure; hold on;
title('All syls, all expts (consol start mean FF)');
% xlabel('acoustic distance');
% ylabel('learning (rel to target)');

Learning_AllExpts=[];
MotifDist_AllExpts=[];
Similar_AllExpts=[];
PreSylSimilar_AllExpts=[];
LearningSEM_all=[];

lt_subplot(2,4,1); hold on; grid on
title('not in target motif')
xlim([-1 1]);
ylim([-1 1]);
ylabel('learning (rel to target)');

lt_subplot(2,4,2); hold on; grid on
title('in target motif')
xlim([-1 1]);
ylim([-1 1]);
ylabel('learning (rel to target)');

lt_subplot(2,4,3:4); hold on; grid on
title('in target motif');
ylim([-1 1]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,5); hold on; grid on
title('not in target motif')
xlim([-1 1]);
ylim([-0.5 0.5]);
ylabel('learning (rel to target)');

lt_subplot(2,4,6); hold on; grid on
title('in target motif')
xlim([-1 1]);
ylim([-0.5 0.5]);
ylabel('learning (rel to target)');

lt_subplot(2,4,7:8); hold on; grid on
title('in target motif');
ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('learning (rel to target)');

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== PLOT ALL SYLS
        syllist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
%         if SeqDepPitch_AcrossBirds.birds{i}.
        
        for j=1:length(syllist);
            syl=syllist{j};
            
            % skip if this is a target
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue;
            end
            
            % ==== collect data
            if norm_all_to_targ==1;
            learning_rel_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            targ_learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
            
            learning_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem;
            learning_sem=learning_sem/abs(targ_learning);
            
            else
                
             learning_rel_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
             % convert to direction of target learning
             learning_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
             learning_rel_targ=learning_rel_targ*learning_dir;
             
            learning_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.sem;
            end
            
            if learning_rel_targ<-3;
                keyboard
            end
            
            Learning_AllExpts=[Learning_AllExpts learning_rel_targ];
            LearningSEM_all=[LearningSEM_all learning_sem];
            
            % - motif distance
            motif_distance=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ;
            MotifDist_AllExpts=[MotifDist_AllExpts motif_distance];
            
            %             % acoustic distance
            %             eucl_dist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
            
            % is similar?
            is_similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            Similar_AllExpts=[Similar_AllExpts is_similar];
            
            % predecing syllable similar?
            preceding_similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            PreSylSimilar_AllExpts=[PreSylSimilar_AllExpts preceding_similar];
            
            % correlations
%             targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
%             
%             try
%             corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(targsyl);
%             catch err
%                 corr_motif=nan;
%             end
%             corr_song=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
            
            
            
            % subplots, for those on same and those on different motif
            if isnan(motif_distance); % then not on same motif
                subplot(2,4,1);
                
                %             ==== plot
                X=-0.5+rand;
                if is_similar==1;
                    lt_plot(X, learning_rel_targ, {'Color','b'});
                else
                    lt_plot(X, learning_rel_targ, {'Color','r'});
                end
                
                % if presyl is similar, outline in black
                if preceding_similar==1;
                    plot(X, learning_rel_targ, 'ok','MarkerSize',7, 'LineWidth', 2.5);
                end
                
                % write text of correlations
%                 if write_corr==1;
%                 text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
%                 end
%                 
                % write syl
                if write_syl==1;
                text(X, learning_rel_targ, [num2str(i) ' ' syl ' ' targsyl])
                end  
                
            else
                % == plot in "in target" subplot
                subplot(2,4,2);
                %     --- plot
                X=-0.5+rand;
                if is_similar==1;
                    lt_plot(X, learning_rel_targ, {'Color','b'});
                else
                    lt_plot(X, learning_rel_targ, {'Color','r'});
                end
                
                % if presyl is similar, outline in black
                if preceding_similar==1;
                    plot(X, learning_rel_targ, 'ok','MarkerSize',7,'LineWidth', 2.5);
                end
                
                                % write syl
                if write_syl==1;
                text(X, learning_rel_targ, [num2str(i) ' ' syl ' ' targsyl])
                end  

                
                
                % === Plot based on motif distance
                subplot(2,4,3:4);
                %  -- plot
                X=motif_distance-0.2+0.4*rand;
                if is_similar==1;
                    lt_plot(X, learning_rel_targ, {'Color','b'});
                else
                    lt_plot(X, learning_rel_targ, {'Color','r'});
                end
                
                % if presyl is similar, outline in black
                if preceding_similar==1;
                    plot(X, learning_rel_targ, 'ok','MarkerSize',7, 'LineWidth', 2.5);
                end
                
                                % write syl
                if write_syl==1;
                text(X, learning_rel_targ, [num2str(i) ' ' syl ' ' targsyl])
                end  

                
            end
        end
    end
end

                if write_syl==1;
                disp('bird syl targsyl')
                end  



%% ==== PLOT ACROSS EXPERIMENTS
% =================================== Different motif
subplot(2,4,5); hold on;
inds=isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

% all
learning_mean=mean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);

% similar
inds2=similar==1;
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);


% diff
inds2=similar==0;
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);



% =================================== Same Motifs (all combined);
subplot(2,4,6); hold on;
inds=~isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

% all
learning_mean=mean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);

% similar
inds2=similar==1;
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);

% diff
inds2=similar==0;
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);



%% ================================= Same Motifs (based on distance
lt_figure; hold on;
if norm_all_to_targ==1;
    ylabel_text='generalization';
else
    ylabel_text='pitch shift (zscore)';
end

inds=~isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
learning_sem=LearningSEM_all(inds);
similar=Similar_AllExpts(inds);
motif_distance=MotifDist_AllExpts(inds);


distances_that_exist=unique(motif_distance);
plotcols=lt_make_plot_colors(length(distances_that_exist),0,0);

hfigs_all=[];

% ===== ALL SYLS
for i=1:length(distances_that_exist);
    distance=distances_that_exist(i);
    
    % -- get those at this distance
    inds2=motif_distance==distance;
    
    Y_learning=learning(inds2);
    Ylearning_sem=learning_sem(inds2);
    
    % ---- RAW VALS
    hfig1=lt_subplot(2,3,1); hold on;
    title('all syls');
    ylabel(ylabel_text);
    xlim([min(distances_that_exist)-1 max(distances_that_exist)+1]);
    xlabel('dist from targ'); set(gca, 'XTick', distances_that_exist)

    if plot_error_bars_on_raw==0;
        lt_plot(distance, Y_learning, {'Color',plotcols{i}});
    else
        lt_plot(distance, Y_learning, {'Color',plotcols{i},'Errors', Ylearning_sem});
        %     errorbar(distance*ones(1,length(Y_learning)), Y_learning, Ylearning_sem, 'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i}, 'MarkerSize',6)
    end
    lt_plot_zeroline;
    line([0 0], ylim, 'Color','k','LineStyle','--');
    
    % -- MEANS
    hfig2=lt_subplot(2,3,4); hold on;
    title('all syls');
    ylabel(ylabel_text);
    xlim([min(distances_that_exist)-1 max(distances_that_exist)+1]);
    xlabel('dist from targ'); set(gca, 'XTick', distances_that_exist)

    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    lt_plot_bar(distance, Ymean, {'Errors',Ysem,'Color','k'});
    hold on;
        line([0 0], ylim, 'Color','k','LineStyle','--');

end
hfigs_all=[hfigs_all hfig1 hfig2];


% ===== SIMILAR
for i=1:length(distances_that_exist);
    
    distance=distances_that_exist(i);
    inds2=motif_distance==distance & similar==1;
        
    if ~any(inds2)
        continue
    end
    
    % run
    Y_learning=learning(inds2);
    Ylearning_sem=learning_sem(inds2);
    
    % ---- RAW VALS
    hfig1=lt_subplot(2,3,2); hold on;
    title('similar syls'); 
    ylabel(ylabel_text);
    xlim([min(distances_that_exist)-1 max(distances_that_exist)+1]);
    xlabel('dist from targ'); set(gca, 'XTick', distances_that_exist)
    
    if plot_error_bars_on_raw==0;
        lt_plot(distance, Y_learning, {'Color',plotcols{i}});
    else
        lt_plot(distance, Y_learning, {'Color',plotcols{i},'Errors', Ylearning_sem});
        %     errorbar(distance*ones(1,length(Y_learning)), Y_learning, Ylearning_sem, 'o','Color',plotcols{i},'MarkerFaceColor',plotcols{i}, 'MarkerSize',6)
    end
    
    lt_plot_zeroline;
        line([0 0], ylim, 'Color','k','LineStyle','--');

    
    % -- MEANS
    hfig2=lt_subplot(2,3,5); hold on;
    title('similar syls');
    ylabel(ylabel_text);
    xlim([min(distances_that_exist)-1 max(distances_that_exist)+1]);
    xlabel('dist from targ'); set(gca, 'XTick', distances_that_exist)

    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    lt_plot_bar(distance, Ymean, {'Errors',Ysem,'Color','b'});
    hold on;
        line([0 0], ylim, 'Color','k','LineStyle','--');

end

hfigs_all=[hfigs_all hfig1 hfig2];

% ===== DIFFERENT
for i=1:length(distances_that_exist);
    
    distance=distances_that_exist(i);
    inds2=motif_distance==distance & similar==0;
        
    if ~any(inds2)
        continue
    end
    
    % run
    Y_learning=learning(inds2);
    Ylearning_sem=learning_sem(inds2);
    
    % ---- RAW VALS
    hfig1=lt_subplot(2,3,3); hold on;
    title('different syls');
    ylabel(ylabel_text);
    xlim([min(distances_that_exist)-1 max(distances_that_exist)+1]);
    xlabel('dist from targ'); set(gca, 'XTick', distances_that_exist)

    if plot_error_bars_on_raw==0;
        lt_plot(distance, Y_learning, {'Color',plotcols{i}});
    else
        lt_plot(distance, Y_learning, {'Color',plotcols{i},'Errors', Ylearning_sem});
    end
    
    lt_plot_zeroline;
        line([0 0], ylim, 'Color','k','LineStyle','--');

    % -- MEANS
    hfig2=lt_subplot(2,3,6); hold on;
    title('different syls');
    ylabel(ylabel_text);
    xlim([min(distances_that_exist)-1 max(distances_that_exist)+1]);
    xlabel('dist from targ'); set(gca, 'XTick', distances_that_exist)

    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    lt_plot_bar(distance, Ymean, {'Errors',Ysem,'Color','b'});
    hold on;
        line([0 0], ylim, 'Color','k','LineStyle','--');

        % -- significantly diff from 0?
        p=signrank(Y_learning);
%         [~, p]=ttest(Y_learning);
        if p<0.1;
            lt_plot_text(distance, Ymean+1.1*Ymean, ['p=' num2str(p,'%3.2g')], 'r');
        end
        
end

hfigs_all=[hfigs_all hfig1 hfig2];

if equal_y_axis==1;
    lt_plot_equalyaxis(hfigs_all);
    
end

linkaxes(hfigs_all, 'x');


%% OLD WAY - all, similar, diff in one plot
% subplot(2,4,7:8); hold on;
% 
% 
% 
% for i=1:length(distances_that_exist);
%     distance=distances_that_exist(i);
%     
%     % -- get those at this distance
%     inds2=motif_distance==distance;
%     
%     Y_learning=learning(inds2);
%     Y_similar=similar(inds2);
%         
%     % -- all
%     Ymean=mean(Y_learning);
%     Ysem=lt_sem(Y_learning);
%     
% %     errorbar(distance-0.1, Ymean, Ysem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);
%     lt_plot_bar(distance-0.1, Ymean, {'Errors',Ysem,'Color','k'});
%     hold on;
%     
%     % -- similar
%     inds3=Y_similar==1;
%     Ymean=mean(Y_learning(inds3));
%     Ysem=lt_sem(Y_learning(inds3));
%     
% %     errorbar(distance, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
%     lt_plot_bar(distance, Ymean, {'Errors',Ysem,'Color','b'});
%     hold on;
%    
%     % -- diff
%     inds3=Y_similar==0;
%     Ymean=mean(Y_learning(inds3));
%     Ysem=lt_sem(Y_learning(inds3));
%     
% %     errorbar(distance+0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
%     lt_plot_bar(distance+0.1, Ymean, {'Errors',Ysem,'Color','r'});
%     hold on;
%     
% end