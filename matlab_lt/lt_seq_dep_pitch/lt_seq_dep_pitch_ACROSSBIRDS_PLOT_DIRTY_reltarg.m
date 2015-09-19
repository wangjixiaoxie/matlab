function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_DIRTY(SeqDepPitch_AcrossBirds, PARAMS)

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% [DIRTY!!!!!!!!!!!!] PLOT CORRELATIONS REL TO DISTANCE IN MOTIF (NOTE: THIS USES THE LEARNING CODE ABOVE, AND IS NOT GREAT. LITERALLY PUTS CORREALTION IN THE LAERNING VARIABLE

write_corr=0;
write_syl=0;
use_motif_corr=0; % else will use song corr

lt_figure; hold on;

Learning_AllExpts=[];
MotifDist_AllExpts=[];
Similar_AllExpts=[];
PreSylSimilar_AllExpts=[];

lt_subplot(2,4,1); hold on; grid on
title('not in target motif')
xlim([-1 1]);
ylim([-1 1]);
ylabel('pair-wise corr (song)');

lt_subplot(2,4,2); hold on; grid on
title('in target motif')
xlim([-1 1]);
ylim([-1 1]);
xlabel('motif distance');
ylabel('pair-wise corr (song)');

lt_subplot(2,4,3:4); hold on; grid on
title('in target motif');
ylim([-1 1]);
xlabel('motif distance');
ylabel('pair-wise corr (song)');

lt_subplot(2,4,5); hold on; grid on
title('not in target motif')
xlim([-1 1]);
ylim([-0.5 0.5]);
ylabel('pair-wise corr (song)');

lt_subplot(2,4,6); hold on; grid on
title('in target motif')
xlim([-1 1]);
ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('pair-wise corr (song)');

lt_subplot(2,4,7:8); hold on; grid on
title('in target motif');
ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('pair-wise corr (song)');

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== PLOT ALL SYLS
        syllist=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed);
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        for j=1:length(syllist);
            syl=syllist{j};
            
            % skip if this is a target
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue;
            end
            
            % ==== collect data
            if use_motif_corr==1;
                try
                    learning_rel_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(targsyl);
                catch err
                    % not in motif
                    learning_rel_targ=nan;
                end
            else
                learning_rel_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
            end
            
            Learning_AllExpts=[Learning_AllExpts learning_rel_targ];
     
%             learning_rel_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.consolid_start_rel_targ;
%             Learning_AllExpts=[Learning_AllExpts learning_rel_targ];
            
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
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            
            try
            corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(targsyl);
            catch err
                corr_motif=nan;
            end
            corr_song=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
            
            
            
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
                if write_corr==1;
                text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                end
                
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
                
                % write text of correlations
                if write_corr==1;
                text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
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
                
                % write text of correlations
                if write_corr==1;
                text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
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



% ==== PLOT ACROSS EXPERIMENTS
% ---- Different motif
subplot(2,4,5); hold on;
inds=isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.SimilarToTarg=similar;

% all
learning_mean=nanmean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);

% similar
inds2=similar==1;
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.similar.LearningRelTarg=Y;

% diff
inds2=similar==0;
Y=learning(inds2);

Ymean=nanmean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.different.LearningRelTarg=Y;



% ---- Same Motifs (all combined);
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



% ---- Same Motifs (based on distance
subplot(2,4,7:8); hold on;
inds=~isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);
motif_distance=MotifDist_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SimilarToTarg=similar;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.DistanceInMotif=motif_distance;


distances_that_exist=unique(motif_distance);
for i=1:length(distances_that_exist);
    distance=distances_that_exist(i);
    
    % -- get those at this distance
    inds2=motif_distance==distance;
    
    Y_learning=learning(inds2);
    Y_similar=similar(inds2);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.distances_that_exist=distances_that_exist;
    
    % -- all
    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    errorbar(distance-0.1, Ymean, Ysem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.all.LearningRelTarg=Y_learning;
    
    % -- similar
    inds3=Y_similar==1;
    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));
    
    errorbar(distance, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.LearningRelTarg=Y_learning(inds3);

    
    % -- diff
    inds3=Y_similar==0;
    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));
    
    errorbar(distance+0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.LearningRelTarg=Y_learning(inds3);
    
    if use_motif_corr==1;
lt_subtitle('THESE ARE CORRELATIONS! (motif) NOT LEARNING')
    else
lt_subtitle('THESE ARE CORRELATIONS! (song) NOT LEARNING')
    end
    
end

%% === [CORRELATIONS] STATISTICAL TESTS ( ON LAERNING VS. MOTIF STUFF) (1. in vs. out of motif) (2. POSITION IN MOTIF)
disp(' ');
disp('------ Statistical tests on effect of motif and position on learning:');
disp(' ');

% ====================== 1) same vs different motif
% ---- All Syls
LearningVals=Learning_AllExpts; % potentially filter similar vs. different
MotifDistanceVals=MotifDist_AllExpts;

% run
inds=isnan(MotifDistanceVals); % different motif

X=LearningVals(inds); % different motif
X=X(~isnan(X));
Y=LearningVals(~inds); % same motif

% test (not paired)
[p, ~] = ranksum(X, Y);

disp(['all syls, motif effect, p = ' num2str(p)]);
figure; hold on; 
lt_plot(1:length(X), X);
lt_plot(1:length(Y), Y, {'Color','r'});
pause;
close all;

% ---- Similar syls
inds=Similar_AllExpts==1;
LearningVals=Learning_AllExpts(inds); % potentially filter similar vs. different
MotifDistanceVals=MotifDist_AllExpts(inds);

% run
inds=isnan(MotifDistanceVals); % different motif

X=LearningVals(inds); % different motif
X=X(~isnan(X));
Y=LearningVals(~inds); % same motif

% test (not paired)
[p, ~] = ranksum(X, Y);

disp(['similar syls, motif effect, p = ' num2str(p)]);
figure; hold on; 
lt_plot(1:length(X), X);
lt_plot(1:length(Y), Y, {'Color','r'});
pause;
close all;

% ---- Diff syls
inds=Similar_AllExpts==0;
LearningVals=Learning_AllExpts(inds); % potentially filter similar vs. different
MotifDistanceVals=MotifDist_AllExpts(inds);

% run
inds=isnan(MotifDistanceVals); % different motif

X=LearningVals(inds); % different motif
X=X(~isnan(X));
Y=LearningVals(~inds); % same motif

% test (not paired)
[p, ~] = ranksum(X, Y);

disp(['different syls, motif effect, p = ' num2str(p)]);
figure; hold on; 
lt_plot(1:length(X), X);
lt_plot(1:length(Y), Y, {'Color','r'});
pause;
close all;


% =============== 2) ANOVA, IGNORE SIMILARITY (WITHIN SAME MOTIF), 
group_anova={};

% first get inds for only those datapoints in same motif
inds=~isnan(MotifDist_AllExpts);

group_anova{1} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% -------------- ONLY USING SIMILAR SYLS
group_anova={};

% first get inds for only those datapoints in same motif
inds=~isnan(MotifDist_AllExpts);

% group_anova{2} = Similar_AllExpts(inds); % similar or different
group_anova{1} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'similarity', 'motifpos'};
GroupNames={'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% =============== ANOVA, IGNORE SIMILARITY (INCLUDING OTHER MOTIF)
group_anova={};

motif_distance_NaNchanged=MotifDist_AllExpts;
motif_distance_NaNchanged(isnan(motif_distance_NaNchanged))=0; % change from nan to 0
group_anova{1} = motif_distance_NaNchanged; % position on motif (or on different motif)
Y = Learning_AllExpts;

GroupNames={'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% =============== 2) ANOVA, SIMILARITY AND MOTIF (anova))
% -- INCLUDING NOT IN MOTIF
% 1) define groups (group cell, each cell is array of IDs for all syls (for
% membership in that group)
group_anova={};
group_anova{1} = Similar_AllExpts; % similar or different
motif_distance_NaNchanged=MotifDist_AllExpts;
motif_distance_NaNchanged(isnan(motif_distance_NaNchanged))=0; % change from nan to 0
group_anova{2} = motif_distance_NaNchanged; % position on motif (or on different motif)

GroupNames={'similarity', 'motifpos'};
Y = Learning_AllExpts;

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% ---- Only looking within motif (linear)
group_anova={};
% first get inds for only those datapoints in same motif
inds=~isnan(MotifDist_AllExpts);

group_anova{1} = Similar_AllExpts(inds); % similar or different
group_anova{2} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'similarity', 'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% ---- Only looking within motif (interaction (simialrity, position))
group_anova={};
% first get inds for only those datapoints in same motif
inds=~isnan(MotifDist_AllExpts);

group_anova{1} = Similar_AllExpts(inds); % similar or different
group_anova{2} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'similarity', 'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames, 'model','interaction');




%% [DIRTY!!!!!!!!!!!!] PLOT ACOUSTIC DISTANCE REL TO DISTANCE IN MOTIF (NOTE: THIS USES THE LEARNING CODE ABOVE, AND IS NOT GREAT. LITERALLY PUTS ACOUSTIC DISTANCE IN THE LAERNING VARIABLE

write_corr=0;
write_syl=1;

lt_figure; hold on;

Learning_AllExpts=[];
MotifDist_AllExpts=[];
Similar_AllExpts=[];
PreSylSimilar_AllExpts=[];

lt_subplot(2,4,1); hold on; grid on
title('not in target motif')
xlim([-1 1]);
% ylim([-1 1]);
ylabel('acoustic distance');

lt_subplot(2,4,2); hold on; grid on
title('in target motif')
xlim([-1 1]);
% ylim([-1 1]);
xlabel('motif distance');
ylabel('acoustic distance');

lt_subplot(2,4,3:4); hold on; grid on
title('in target motif');
% ylim([-1 1]);
xlabel('motif distance');
ylabel('acoustic distance');

lt_subplot(2,4,5); hold on; grid on
title('not in target motif')
xlim([-1 1]);
% ylim([-0.5 0.5]);
ylabel('acoustic distance');

lt_subplot(2,4,6); hold on; grid on
title('in target motif')
xlim([-1 1]);
% ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('acoustic distance');

lt_subplot(2,4,7:8); hold on; grid on
title('in target motif');
% ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('acoustic distance');

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== PLOT ALL SYLS
        syllist=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed);
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        for j=1:length(syllist);
            syl=syllist{j};
            
            % skip if this is a target
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue;
            end
            
            % ==== collect data
            learning_rel_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
            Learning_AllExpts=[Learning_AllExpts learning_rel_targ];
     
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
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            
            try
            corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(targsyl);
            catch err
                corr_motif=nan;
            end
            corr_song=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
            
            
            
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
                if write_corr==1;
                text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                end
                
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
                
                % write text of correlations
                if write_corr==1;
                text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
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
                
                % write text of correlations
                if write_corr==1;
                text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
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



% ==== PLOT ACROSS EXPERIMENTS
% ---- Different motif
subplot(2,4,5); hold on;
inds=isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.SimilarToTarg=similar;

% all
learning_mean=nanmean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);

% similar
inds2=similar==1;
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.similar.LearningRelTarg=Y;

% diff
inds2=similar==0;
Y=learning(inds2);

Ymean=nanmean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.different.LearningRelTarg=Y;



% ---- Same Motifs (all combined);
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



% ---- Same Motifs (based on distance
subplot(2,4,7:8); hold on;
inds=~isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);
motif_distance=MotifDist_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SimilarToTarg=similar;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.DistanceInMotif=motif_distance;


distances_that_exist=unique(motif_distance);
for i=1:length(distances_that_exist);
    distance=distances_that_exist(i);
    
    % -- get those at this distance
    inds2=motif_distance==distance;
    
    Y_learning=learning(inds2);
    Y_similar=similar(inds2);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.distances_that_exist=distances_that_exist;
    
    % -- all
    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    errorbar(distance-0.1, Ymean, Ysem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.all.LearningRelTarg=Y_learning;
    
    % -- similar
    inds3=Y_similar==1;
    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));
    
    errorbar(distance, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.LearningRelTarg=Y_learning(inds3);

    
    % -- diff
    inds3=Y_similar==0;
    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));
    
    errorbar(distance+0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.LearningRelTarg=Y_learning(inds3);
    
lt_subtitle('THESE ARE ACOUSTIC DISTNACE! NOT LEARNING')
    
end


%% === [ACOUSTIC] STATISTICAL TESTS ( ON LAERNING VS. MOTIF STUFF) (1. in vs. out of motif) (2. POSITION IN MOTIF)
disp(' ');
disp('------ Statistical tests on effect of motif and position on learning:');
disp(' ');

% ====================== 1) same vs different motif
% ---- All Syls
LearningVals=Learning_AllExpts; % potentially filter similar vs. different
MotifDistanceVals=MotifDist_AllExpts;

% run
inds=isnan(MotifDistanceVals); % different motif

X=LearningVals(inds); % different motif
X=X(~isnan(X));
Y=LearningVals(~inds); % same motif

% test (not paired)
[p, ~] = ranksum(X, Y);

disp(['all syls, motif effect, p = ' num2str(p)]);
figure; hold on; 
lt_plot(1:length(X), X);
lt_plot(1:length(Y), Y, {'Color','r'});
pause;
close all;

% ---- Similar syls
inds=Similar_AllExpts==1;
LearningVals=Learning_AllExpts(inds); % potentially filter similar vs. different
MotifDistanceVals=MotifDist_AllExpts(inds);

% run
inds=isnan(MotifDistanceVals); % different motif

X=LearningVals(inds); % different motif
X=X(~isnan(X));
Y=LearningVals(~inds); % same motif

% test (not paired)
[p, ~] = ranksum(X, Y);

disp(['similar syls, motif effect, p = ' num2str(p)]);
figure; hold on; 
lt_plot(1:length(X), X);
lt_plot(1:length(Y), Y, {'Color','r'});
pause;
close all;

% ---- Diff syls
inds=Similar_AllExpts==0;
LearningVals=Learning_AllExpts(inds); % potentially filter similar vs. different
MotifDistanceVals=MotifDist_AllExpts(inds);

% run
inds=isnan(MotifDistanceVals); % different motif

X=LearningVals(inds); % different motif
X=X(~isnan(X));
Y=LearningVals(~inds); % same motif

% test (not paired)
[p, ~] = ranksum(X, Y);

disp(['different syls, motif effect, p = ' num2str(p)]);
figure; hold on; 
lt_plot(1:length(X), X);
lt_plot(1:length(Y), Y, {'Color','r'});
pause;
close all;


% =============== 2) ANOVA, IGNORE SIMILARITY (WITHIN SAME MOTIF), 
group_anova={};

% first get inds for only those datapoints in same motif
inds=~isnan(MotifDist_AllExpts);

group_anova{1} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% -------------- ONLY USING SIMILAR SYLS
group_anova={};

% first get inds for only those datapoints in same motif
inds=~isnan(MotifDist_AllExpts);

% group_anova{2} = Similar_AllExpts(inds); % similar or different
group_anova{1} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'similarity', 'motifpos'};
GroupNames={'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% =============== ANOVA, IGNORE SIMILARITY (INCLUDING OTHER MOTIF)
group_anova={};

motif_distance_NaNchanged=MotifDist_AllExpts;
motif_distance_NaNchanged(isnan(motif_distance_NaNchanged))=0; % change from nan to 0
group_anova{1} = motif_distance_NaNchanged; % position on motif (or on different motif)
Y = Learning_AllExpts;

GroupNames={'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% =============== 2) ANOVA, SIMILARITY AND MOTIF (anova))
% -- INCLUDING NOT IN MOTIF
% 1) define groups (group cell, each cell is array of IDs for all syls (for
% membership in that group)
group_anova={};
group_anova{1} = Similar_AllExpts; % similar or different
motif_distance_NaNchanged=MotifDist_AllExpts;
motif_distance_NaNchanged(isnan(motif_distance_NaNchanged))=0; % change from nan to 0
group_anova{2} = motif_distance_NaNchanged; % position on motif (or on different motif)

GroupNames={'similarity', 'motifpos'};
Y = Learning_AllExpts;

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% ---- Only looking within motif (linear)
group_anova={};
% first get inds for only those datapoints in same motif
inds=~isnan(MotifDist_AllExpts);

group_anova{1} = Similar_AllExpts(inds); % similar or different
group_anova{2} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'similarity', 'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames);


% ---- Only looking within motif (interaction (simialrity, position))
group_anova={};
% first get inds for only those datapoints in same motif
inds=~isnan(MotifDist_AllExpts);

group_anova{1} = Similar_AllExpts(inds); % similar or different
group_anova{2} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'similarity', 'motifpos'};

[p, table, stats, terms]=anovan(Y, group_anova, 'varnames', GroupNames, 'model','interaction');



%% %% [DIRTY!!!!!!!!!!!!] [PAIRWISE STATISTICS] [PLOT ACOUSTIC DISTANCE REL TO DISTANCE IN MOTIF (NOTE: THIS USES THE LEARNING CODE ABOVE, AND IS NOT GREAT. LITERALLY PUTS ACOUSTIC DISTANCE IN THE LAERNING VARIABLE
clear syl;
write_corr=0;
write_syl=1;

% ONLY ONE OF THESE CAN BE 1
plot_acoustic=1;
plot_corr_song=0;
plot_corr_motif=0;


lt_figure; hold on;

Learning_AllExpts=[];
MotifDist_AllExpts=[];
Similar_AllExpts=[];
PreSylSimilar_AllExpts=[];

lt_subplot(2,4,1); hold on; grid on
title('not in target motif')
xlim([-1 1]);
% ylim([-1 1]);
ylabel('learning (rel to target)');

lt_subplot(2,4,2); hold on; grid on
title('in target motif')
xlim([-1 1]);
% ylim([-1 1]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,3:4); hold on; grid on
title('in target motif');
% ylim([-1 1]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,5); hold on; grid on
title('not in target motif')
xlim([-1 1]);
% ylim([-0.5 0.5]);
ylabel('learning (rel to target)');

lt_subplot(2,4,6); hold on; grid on
title('in target motif')
xlim([-1 1]);
% ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,7:8); hold on; grid on
title('in target motif');
% ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('learning (rel to target)');

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== PLOT ALL SYLS
        syllist=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed);
        %         targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        for j=1:length(syllist);
            syl1=syllist{j};
            % skip is last syl
            if j==length(syllist);
                continue;
            end
            
            for jj=j+1:length(syllist);
                syl2=syllist{jj};
                
                % skip if this is a target
                %             if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                %                 continue;
                %             end
                
                % ==== collect data
                % acoustic distance
                fv1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl1).fv_baseline_zscore_mean;
                fv2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl2).fv_baseline_zscore_mean;
                acoustdist=sqrt(sum((fv1-fv2).^2));

                % correlations
                try
                    corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                catch err
                    corr_motif=nan;
                end
                corr_song=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2);
                
                
                % put into variable that will plot
                if plot_acoustic==1;
                    learning_rel_targ=acoustdist;
                elseif plot_corr_song==1;
                    learning_rel_targ=corr_song;
                elseif plot_corr_motif==1;
                    learning_rel_targ=corr_motif;
                end
                Learning_AllExpts=[Learning_AllExpts learning_rel_targ];
                
                
                % - motif distance
                motifnum1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).regexp_motifnum;
                motifnum2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_motifnum;
                if motifnum1~=motifnum2;
                    motif_distance=nan;
                else
                    motifpos1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).regexp_PosInMotif_thissyl;
                    motifpos2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_PosInMotif_thissyl;
                    
                    motif_distance=abs(motifpos1-motifpos2);
                end
                
                MotifDist_AllExpts=[MotifDist_AllExpts motif_distance];
                
                %             % acoustic distance
                %             eucl_dist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
                
                % is similar?
                lowersyl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).single_syl;
                lowersyl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).single_syl;
                if lowersyl1==lowersyl2;
                    is_similar=1;
                else
                    is_similar=0;
                end
                
                Similar_AllExpts=[Similar_AllExpts is_similar];
                
                
                % predecing syllable similar?
                presyl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).preceding_syl;
                if length(presyl1)>1;
                    presyl1=regexp(presyl1, '[A-Z]', 'match');
                    presyl1=lower(presyl1);
                end
                presyl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).preceding_syl;
                if length(presyl2)>1;
                    presyl2=regexp(presyl2, '[A-Z]', 'match');
                    presyl2=lower(presyl2);
                end
                
                if strcmp(presyl1, presyl2);
                    preceding_similar=1;
                else
                    preceding_similar=0;
                end
                PreSylSimilar_AllExpts=[PreSylSimilar_AllExpts preceding_similar];
                
                
                
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
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
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
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
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
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
                    end
                    
                    
                end
            end
        end
    end
end

if write_syl==1;
    disp('bird syl1 syl2')
end



% ==== PLOT ACROSS EXPERIMENTS
% ---- Different motif
subplot(2,4,5); hold on;
inds=isnan(MotifDist_AllExpts);

preceding=PreSylSimilar_AllExpts(inds);
learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.SimilarToTarg=similar;

% all
learning_mean=nanmean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);


% similar
inds2=similar==1;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.similar.LearningRelTarg=Y;

% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerSize', 8);


% diff
inds2=similar==0;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=nanmean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.different.LearningRelTarg=Y;
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=nanmean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerSize', 8);




% ---- Same Motifs (all combined);
subplot(2,4,6); hold on;
inds=~isnan(MotifDist_AllExpts);

preceding=PreSylSimilar_AllExpts(inds);
learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

% all
learning_mean=mean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);

% similar
inds2=similar==1;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerSize', 8);

% diff
inds2=similar==0;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerSize', 8);



% ---- Same Motifs (based on distance
subplot(2,4,7:8); hold on;
inds=~isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);
motif_distance=MotifDist_AllExpts(inds);
preceding=PreSylSimilar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SimilarToTarg=similar;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.DistanceInMotif=motif_distance;


distances_that_exist=unique(motif_distance);
for i=1:length(distances_that_exist);
    distance=distances_that_exist(i);
    
    % -- get those at this distance
    inds2=motif_distance==distance;
    Y_preceding=preceding(inds2);
    Y_learning=learning(inds2);
    Y_similar=similar(inds2);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.distances_that_exist=distances_that_exist;
    
    % -- all
    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    errorbar(distance-0.1, Ymean, Ysem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.all.LearningRelTarg=Y_learning;
    
    % -- similar
    inds3=Y_similar==1;
    Y_pre=Y_preceding(inds3);

    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));
    
    errorbar(distance, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.LearningRelTarg=Y_learning(inds3);
    
    Y=Y_learning(inds3);
    % plot means with presyl similar/diff
    inds4=Y_pre==1;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
    % plot means with presyl similar/diff
    inds4=Y_pre==0;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'ob', 'MarkerSize', 8);
    
    
    % -- diff
    inds3=Y_similar==0;
    Y_pre=Y_preceding(inds3);
    
    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));

    errorbar(distance+0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.LearningRelTarg=Y_learning(inds3);
        
    Y=Y_learning(inds3);
    % plot means with presyl similar/diff
    inds4=Y_pre==1;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
    % plot means with presyl similar/diff
    inds4=Y_pre==0;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'or', 'MarkerSize', 8);

    
    
    
    if plot_acoustic==1;
        lt_subtitle('ALL PAIRS OF SYLS - ACOUSTIC vs. MOTIF POSITION')
    elseif plot_corr_song==1;
        lt_subtitle('ALL PAIRS OF SYLS - CORR (song) vs. MOTIF POSITION')
    elseif plot_corr_motif==1;
        lt_subtitle('ALL PAIRS OF SYLS - COPR (motif) vs. MOTIF POSITION')
    end
    
end

%% %% [DIRTY!!!!!!!!!!!!] [PAIRWISE STATISTICS] [PLOT CORR DISTANCE REL TO DISTANCE IN MOTIF (NOTE: THIS USES THE LEARNING CODE ABOVE, AND IS NOT GREAT. LITERALLY PUTS ACOUSTIC DISTANCE IN THE LAERNING VARIABLE
clear syl;
write_corr=0;
write_syl=1;

% ONLY ONE OF THESE CAN BE 1
plot_acoustic=0;
plot_corr_song=1;
plot_corr_motif=0;


lt_figure; hold on;

Learning_AllExpts=[];
MotifDist_AllExpts=[];
Similar_AllExpts=[];
PreSylSimilar_AllExpts=[];

lt_subplot(2,4,1); hold on; grid on
title('not in target motif')
xlim([-1 1]);
% ylim([-1 1]);
ylabel('learning (rel to target)');

lt_subplot(2,4,2); hold on; grid on
title('in target motif')
xlim([-1 1]);
% ylim([-1 1]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,3:4); hold on; grid on
title('in target motif');
% ylim([-1 1]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,5); hold on; grid on
title('not in target motif')
xlim([-1 1]);
% ylim([-0.5 0.5]);
ylabel('learning (rel to target)');

lt_subplot(2,4,6); hold on; grid on
title('in target motif')
xlim([-1 1]);
% ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,7:8); hold on; grid on
title('in target motif');
% ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('learning (rel to target)');

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== PLOT ALL SYLS
        syllist=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed);
        %         targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        for j=1:length(syllist);
            syl1=syllist{j};
            % skip is last syl
            if j==length(syllist);
                continue;
            end
            
            for jj=j+1:length(syllist);
                syl2=syllist{jj};
                
                % skip if this is a target
                %             if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                %                 continue;
                %             end
                
                % ==== collect data
                % acoustic distance
                fv1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl1).fv_baseline_zscore_mean;
                fv2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl2).fv_baseline_zscore_mean;
                acoustdist=sqrt(sum((fv1-fv2).^2));

                % correlations
                try
                    corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl2);
                catch err
                    corr_motif=nan;
                end
                corr_song=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).CORRELATIONS.song_by_song.corrcoeff_vs.(syl2);
                
                
                % put into variable that will plot
                if plot_acoustic==1;
                    learning_rel_targ=acoustdist;
                elseif plot_corr_song==1;
                    learning_rel_targ=corr_song;
                elseif plot_corr_motif==1;
                    learning_rel_targ=corr_motif;
                end
                Learning_AllExpts=[Learning_AllExpts learning_rel_targ];
                
                
                % - motif distance
                motifnum1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).regexp_motifnum;
                motifnum2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_motifnum;
                if motifnum1~=motifnum2;
                    motif_distance=nan;
                else
                    motifpos1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).regexp_PosInMotif_thissyl;
                    motifpos2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).regexp_PosInMotif_thissyl;
                    
                    motif_distance=abs(motifpos1-motifpos2);
                end
                
                MotifDist_AllExpts=[MotifDist_AllExpts motif_distance];
                
                %             % acoustic distance
                %             eucl_dist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
                
                % is similar?
                lowersyl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).single_syl;
                lowersyl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).single_syl;
                if lowersyl1==lowersyl2;
                    is_similar=1;
                else
                    is_similar=0;
                end
                
                Similar_AllExpts=[Similar_AllExpts is_similar];
                
                
                % predecing syllable similar?
                presyl1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl1).preceding_syl;
                if length(presyl1)>1;
                    presyl1=regexp(presyl1, '[A-Z]', 'match');
                    presyl1=lower(presyl1);
                end
                presyl2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl2).preceding_syl;
                if length(presyl2)>1;
                    presyl2=regexp(presyl2, '[A-Z]', 'match');
                    presyl2=lower(presyl2);
                end
                
                if strcmp(presyl1, presyl2);
                    preceding_similar=1;
                else
                    preceding_similar=0;
                end
                PreSylSimilar_AllExpts=[PreSylSimilar_AllExpts preceding_similar];
                
                
                
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
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
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
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
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
                    
                    % write text of correlations
                    if write_corr==1;
                        text(X, learning_rel_targ, [num2str(corr_motif, '%2.1f') ' ' num2str(corr_song, '%2.1f')])
                    end
                    
                    % write syl
                    if write_syl==1;
                        text(X, learning_rel_targ, [num2str(i) ' ' syl1 ' ' syl2])
                    end
                    
                    
                end
            end
        end
    end
end

if write_syl==1;
    disp('bird syl1 syl2')
end



% ==== PLOT ACROSS EXPERIMENTS
% ---- Different motif
subplot(2,4,5); hold on;
inds=isnan(MotifDist_AllExpts);

preceding=PreSylSimilar_AllExpts(inds);
learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.all.SimilarToTarg=similar;

% all
learning_mean=nanmean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);


% similar
inds2=similar==1;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.similar.LearningRelTarg=Y;

% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerSize', 8);


% diff
inds2=similar==0;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=nanmean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.diff_motif.different.LearningRelTarg=Y;
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=nanmean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerSize', 8);




% ---- Same Motifs (all combined);
subplot(2,4,6); hold on;
inds=~isnan(MotifDist_AllExpts);

preceding=PreSylSimilar_AllExpts(inds);
learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);

% all
learning_mean=mean(learning);
learning_sem=lt_sem(learning);

errorbar(-0.1, learning_mean, learning_sem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);

% similar
inds2=similar==1;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'ob', 'MarkerSize', 8);

% diff
inds2=similar==0;
Y_pre=preceding(inds2);
Y=learning(inds2);

Ymean=mean(Y);
Ysem=lt_sem(Y);
errorbar(0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
% plot means with presyl similar/diff
inds4=Y_pre==1;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
% plot means with presyl similar/diff
inds4=Y_pre==0;
Ymean=mean(Y(inds4));
Ysem=lt_sem(Y(inds4));
errorbar(0+0.5, Ymean, Ysem, 'or', 'MarkerSize', 8);



% ---- Same Motifs (based on distance
subplot(2,4,7:8); hold on;
inds=~isnan(MotifDist_AllExpts);

learning=Learning_AllExpts(inds);
similar=Similar_AllExpts(inds);
motif_distance=MotifDist_AllExpts(inds);
preceding=PreSylSimilar_AllExpts(inds);

SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.LearningRelTarg=learning;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SimilarToTarg=similar;
SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.DistanceInMotif=motif_distance;


distances_that_exist=unique(motif_distance);
for i=1:length(distances_that_exist);
    distance=distances_that_exist(i);
    
    % -- get those at this distance
    inds2=motif_distance==distance;
    Y_preceding=preceding(inds2);
    Y_learning=learning(inds2);
    Y_similar=similar(inds2);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.distances_that_exist=distances_that_exist;
    
    % -- all
    Ymean=mean(Y_learning);
    Ysem=lt_sem(Y_learning);
    
    errorbar(distance-0.1, Ymean, Ysem, 'sk', 'MarkerFaceColor','k', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.all.LearningRelTarg=Y_learning;
    
    % -- similar
    inds3=Y_similar==1;
    Y_pre=Y_preceding(inds3);

    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));
    
    errorbar(distance, Ymean, Ysem, 'sb', 'MarkerFaceColor','b', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.similar.LearningRelTarg=Y_learning(inds3);
    
    Y=Y_learning(inds3);
    % plot means with presyl similar/diff
    inds4=Y_pre==1;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'ob', 'MarkerFaceColor','k', 'MarkerSize', 8);
    % plot means with presyl similar/diff
    inds4=Y_pre==0;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'ob', 'MarkerSize', 8);
    
    
    % -- diff
    inds3=Y_similar==0;
    Y_pre=Y_preceding(inds3);
    
    Ymean=mean(Y_learning(inds3));
    Ysem=lt_sem(Y_learning(inds3));

    errorbar(distance+0.1, Ymean, Ysem, 'sr', 'MarkerFaceColor','r', 'MarkerSize', 10);
    
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.distance=distance;
    SeqDepPitch_AcrossBirds.SUMMARY_DATA.MOTIFS.data.same_motif.SORTED_BY_DISTANCE{i}.different.LearningRelTarg=Y_learning(inds3);
        
    Y=Y_learning(inds3);
    % plot means with presyl similar/diff
    inds4=Y_pre==1;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'or', 'MarkerFaceColor','k', 'MarkerSize', 8);
    % plot means with presyl similar/diff
    inds4=Y_pre==0;
    Ymean=mean(Y(inds4));
    Ysem=lt_sem(Y(inds4));
    errorbar(distance+0.3, Ymean, Ysem, 'or', 'MarkerSize', 8);

    
    
    
    if plot_acoustic==1;
        lt_subtitle('ALL PAIRS OF SYLS - ACOUSTIC vs. MOTIF POSITION')
    elseif plot_corr_song==1;
        lt_subtitle('ALL PAIRS OF SYLS - CORR (song) vs. MOTIF POSITION')
    elseif plot_corr_motif==1;
        lt_subtitle('ALL PAIRS OF SYLS - COPR (motif) vs. MOTIF POSITION')
    end
    
end

