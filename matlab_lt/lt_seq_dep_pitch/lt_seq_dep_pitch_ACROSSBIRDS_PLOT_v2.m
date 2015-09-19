function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_PLOT_v2(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 8/13/15 - v2, uses learning metric, not consol start FF.

%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% PLOT CORRELATION POPULATION INFORMATION

% 1) For each experiment, one plot, for each syl, plot bar plot against all
% other syls. (song-by-song correlation)
for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        SylFields_Unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        lt_figure; hold on;
        
        for j=1:length(SylFields_Unique);
            syl=SylFields_Unique{j};
            
            lt_subplot(ceil(length(SylFields_Unique)/3), 3, j); hold on; title(syl);
            
            Y_corr=[];
            X=[];
            
            % ---- GET CORR AGAINST OTHER SYLS (SONG BY SONG)
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                X=[X jj];
                Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(othersyl)];
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(othersyl));
                    Y_corr=[Y_corr nan];
                end
            end
            
            % PLOT
            bar(X-0.2, Y_corr, 0.5);
            %             xlim([0.5 length(SylFields_Unique)+0.5]);
            ylim([-0.2, 1]);
            set(gca, 'XTick', 1:length(SylFields_Unique))
            set(gca, 'XTickLabel', SylFields_Unique);
            
            % === OVERLAY OLD DATA FORMAT, TO MAKE SURE EXTRACTED DATA ARE
            % CORRECT
            ind_of_this_syl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.ind_of_this_syl;
            Y_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.RhoMat(:,ind_of_this_syl);
            %             Y_old(isnan(Y_old))=[];
            lt_plot(1:length(Y_old), Y_old, {'Color','r'});
            
            
            % ============== GET CORR WITHIN MOTIF
            X=[];
            Y_corr=[];
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                
                % make sure this syl is in the same motif
                
                X=[X jj];
                try
                    Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(othersyl)];
                catch err
                    Y_corr=[Y_corr nan];
                end
                
            end
            
            bar(X+0.2, Y_corr, 0.5, 'FaceColor', [0.1 0.9 0.1]);
            
            
            % === OVERLAY OLD DATA FORMAT, TO MAKE SURE EXTRACTED DATA ARE
            % CORRECT
            try
                motif_thissyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.regexp_motifnum;
                pos_in_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.regexp_PosInMotif;
                
                Y_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{motif_thissyl}.sub_class{1}.CORRELATIONS.SUBTRACT_SONG_MEAN.RhoMat(pos_in_motif,:);
                
                %             Y_old(isnan(Y_old))=[];
                lt_plot(1:length(Y_old), Y_old, {'Color','g'});
            catch err
            end
            
        end
        
        lt_subtitle(['Correlations (bars=extracted, dots=old(order might be wrong); Bird: ' birdname ', expt ' exptname]);
    end
end
pause;
close all;



%% 2) PLOT PAIRWISE CORRELATION STATS (across all expts and birds)
AllExptsData=struct;
AllExptsData.pairwise_corr=[];
AllExptsData.pairwise_corr_pval=[];
AllExptsData.pair_is_similar=[];

AllExptsData.AcousticDist_all=[];

AllExptsData.CorrCoeff_motif_all=[];
AllExptsData.Pval_motif_all=[];

AllExptsData.CorrCoeff_motif_SubtrSong_all=[];
AllExptsData.Pval_motif_SubtrSong_all=[];

for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    hfig1=lt_figure; hold on;
    hfig2=lt_figure; hold on;
    hfig3=lt_figure; hold on;
    
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        SylFields_Unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        CorrCoeff_all=[];
        Pval_all=[];
        Similar_status_all=[];
        AcousticDist_all=[];
        
        CorrCoeff_motif_all=[];
        Pval_motif_all=[];
        
        CorrCoeff_motif_SubtrSong_all=[];
        Pval_motif_SubtrSong_all=[];
        
        for j=1:length(SylFields_Unique);
            syl=SylFields_Unique{j};
            
            % =========  GET CORR AGAINST ALL SYLLABLES THAT COME AFTER IT (no double
            % counting, and no auto-corr).
            if j==length(SylFields_Unique); % last syl - stop.
                continue;
            end
            
            for jj=j+1:length(SylFields_Unique);
                
                othersyl=SylFields_Unique{jj};
                
                
                % ====== EXTRACT DATA
                % -- Corr (song-by_song)
                corr=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(othersyl);
                p_val=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.p_val_vs.(othersyl);
                
                % --- Corr (motif_by_motif)
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs, othersyl); % might not be in same motif
                    corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(othersyl);
                    pval_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.p_val_vs.(othersyl);
                    
                    corr_motif_MinusSongMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(othersyl);
                    pval_motif_MinusSongMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.p_val_vs.(othersyl);
                else
                    corr_motif=nan;
                    pval_motif=nan;
                    
                    corr_motif_MinusSongMean=nan;
                    pval_motif_MinusSongMean=nan;
                end
                
                % -- Similar vs. different
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl == SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).single_syl
                    issimilar=1;
                else
                    issimilar=0;
                end
                
                % --- Motif distance (nan = different) - in progress
                
                
                
                % --- Get acoustic difference - IN PROGRESS
                fv1_zscore=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_zscore_mean;
                fv2_zscore=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(othersyl).fv_baseline_zscore_mean;
                acoustic_dist=sqrt(sum((fv1_zscore-fv2_zscore).^2));
                
                % --- OUTPUT
                if ~isempty(corr) && ~isnan(corr);
                    CorrCoeff_all=[CorrCoeff_all corr];
                    Pval_all=[Pval_all p_val];
                    Similar_status_all=[Similar_status_all issimilar];
                    AcousticDist_all=[AcousticDist_all acoustic_dist];
                    
                    CorrCoeff_motif_all=[CorrCoeff_motif_all corr_motif];
                    Pval_motif_all=[Pval_motif_all pval_motif];
                    
                    CorrCoeff_motif_SubtrSong_all=[CorrCoeff_motif_SubtrSong_all corr_motif_MinusSongMean];
                    Pval_motif_SubtrSong_all=[Pval_motif_SubtrSong_all pval_motif_MinusSongMean];
                    
                    
                    % ===== OUTPUT FOR ALL DATA ================
                    AllExptsData.pairwise_corr=[AllExptsData.pairwise_corr corr];
                    AllExptsData.pairwise_corr_pval=[AllExptsData.pairwise_corr_pval p_val];
                    AllExptsData.pair_is_similar=[AllExptsData.pair_is_similar issimilar];
                    
                    AllExptsData.AcousticDist_all=[AllExptsData.AcousticDist_all acoustic_dist];
                    
                    AllExptsData.CorrCoeff_motif_all=[AllExptsData.CorrCoeff_motif_all corr_motif];
                    AllExptsData.Pval_motif_all=[AllExptsData.Pval_motif_all pval_motif];
                    
                    AllExptsData.CorrCoeff_motif_SubtrSong_all=[AllExptsData.CorrCoeff_motif_SubtrSong_all corr_motif_MinusSongMean];
                    AllExptsData.Pval_motif_SubtrSong_all=[AllExptsData.Pval_motif_SubtrSong_all pval_motif_MinusSongMean];
                    
                    % =========================================
                end
                
            end
            
        end
        
        % ==== PLOT (for this bird/epxt);
        figure(hfig1);
        
        lt_subplot(ceil(NumExperiments/2), 2, ii); % for this experiment
        title(exptname); hold on;
        
        % 1a) pval vs. corr (all, similar, diff) (song-by-song
        xlabel('corr coeff');
        ylabel('log10 p-value of correlation');
        grid on;
        
        % similar
        inds=Similar_status_all==1;
        X=CorrCoeff_all(inds);
        Y=Pval_all(inds);
        plot(X, log10(Y), 'ob');
        plot(mean(X), log10(mean(Y)), 'ob', 'MarkerSize', 12)
        
        % diff
        inds=Similar_status_all==0;
        X=CorrCoeff_all(inds);
        Y=Pval_all(inds);
        plot(X, log10(Y), 'or');
        plot(mean(X), log10(mean(Y)), 'or', 'MarkerSize', 12)
        
        % 1b) same, but motif by motif
        % similar
        inds=Similar_status_all==1;
        X=CorrCoeff_motif_all(inds);
        Y=Pval_motif_all(inds);
        plot(X, log10(Y), 'sb');
        
        plot(nanmean(X), log10(nanmean(Y)), 'sb', 'MarkerSize', 12)
        
        % diff
        inds=Similar_status_all==0;
        X=CorrCoeff_motif_all(inds);
        Y=Pval_motif_all(inds);
        plot(X, log10(Y), 'sr');
        
        plot(nanmean(X), log10(nanmean(Y)), 'sr', 'MarkerSize', 12)
        
        % 1c) same, but motif by motif (song-mean subtracted)
        % similar
        inds=Similar_status_all==1;
        X=CorrCoeff_motif_SubtrSong_all(inds);
        Y=Pval_motif_SubtrSong_all(inds);
        plot(X, log10(Y), '.b');
        
        plot(nanmean(X), log10(nanmean(Y)), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 8)
        
        % diff
        inds=Similar_status_all==0;
        X=CorrCoeff_motif_SubtrSong_all(inds);
        Y=Pval_motif_SubtrSong_all(inds);
        plot(X, log10(Y), '.r');
        
        plot(nanmean(X), log10(nanmean(Y)), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
        
        
        % OTHER THINGS
        line(xlim, [log10(0.05) log10(0.05)], 'Color', 'k')
        
        
        % 2) Do the different methods of calcualting correaltions
        % correlate?
        figure(hfig2);
        
        lt_subplot(ceil(NumExperiments/2), 2, ii); % for this experiment
        title(exptname); hold on;
        grid on; box on
        
        % 2a) song-by-song vs. motif-by-motif
        xlabel('motif-by-motif corr');
        ylabel('song-by-song corr');
        
        % similar
        inds=Similar_status_all==1;
        X=CorrCoeff_motif_all(inds);
        Y=CorrCoeff_all(inds);
        
        inds2=~isnan(X); % remove nans
        X=X(inds2);
        Y=Y(inds2);
        
        lt_plot(X,Y, {'Color','b'});
        
        % different
        inds=Similar_status_all==0;
        X=CorrCoeff_motif_all(inds);
        Y=CorrCoeff_all(inds);
        
        inds2=~isnan(X); % remove nans
        X=X(inds2);
        Y=Y(inds2);
        
        lt_plot(X,Y, {'Color','r'});
        
        
        % 2b) song-by-song vs. motif-by-motif(minus song)
        xlabel('motif-by-motif corr (open = song mean subtracted)');
        ylabel('song-by-song corr');
        
        % similar
        inds=Similar_status_all==1;
        X=CorrCoeff_motif_SubtrSong_all(inds);
        Y=CorrCoeff_all(inds);
        
        inds2=~isnan(X); % remove nans
        X=X(inds2);
        Y=Y(inds2);
        
        plot(X,Y, 'ob');
        
        % diff
        inds=Similar_status_all==0;
        X=CorrCoeff_motif_SubtrSong_all(inds);
        Y=CorrCoeff_all(inds);
        
        inds2=~isnan(X); % remove nans
        X=X(inds2);
        Y=Y(inds2);
        
        plot(X,Y, 'or');
        
        % Other stuff
        xlim([-0.3 1]);
        ylim([-0.3 1]);
        
        line([0 0], ylim,'Color', 'k');
        line(xlim, [0 0],'Color', 'k');
        line(xlim, ylim, 'Color','k')
        
        
        % 3) Do pairwise correlations correlate with acoustic distance?
        figure(hfig3); hold on;
        lt_subplot(ceil(NumExperiments/2), 2, ii); % for this experiment
        title([exptname ' (regression is using similar only)']); hold on;
        grid on; box on
        
        % 3a) acoustic dist vs. song by song corr
        xlabel('corr (closed song by song) (open=motifbymotif, song subtr)');
        ylabel('acoustic euclidian dist');
        
        % similar
        inds=Similar_status_all==1;
        X=CorrCoeff_all(inds);
        Y=AcousticDist_all(inds);
        
        lt_regress(Y,X, 1);
        
        % different
        inds=Similar_status_all==0;
        X=CorrCoeff_all(inds);
        Y=AcousticDist_all(inds);
        
        lt_plot(X, Y, {'Color', 'r'});
        
        % 3b) same, but plotting motifbymotif corr
        % similar
        inds=Similar_status_all==1;
        X=CorrCoeff_motif_SubtrSong_all(inds);
        Y=AcousticDist_all(inds);
        
        plot(X, Y, 'ob');
        
        % different
        inds=Similar_status_all==0;
        X=CorrCoeff_motif_SubtrSong_all(inds);
        Y=AcousticDist_all(inds);
        
        plot(X, Y, 'or');
        
        % 4) Motif position/(same or diff motif) - IN PROGRESS
        
    end
    
    
    figure(hfig1);
    lt_subtitle(['pairwise correlations (circle=song; squ=motif, dot=motif, minus song) - ' birdname])
    
    figure(hfig2);
    lt_subtitle(['pairwise correlations, comparing different ways of calcualting - ' birdname])
    
    figure(hfig3);
    lt_subtitle(['acoustic distance vs. corr- ' birdname])
    
end
pause; close all;


%% THINGS TO PLOT (across all PAIRS)
% 1) motif vs. song corr similar?
% 2) acoustic distance vs corr?
% 3) motif position/same motif, diff?
% 4) similarity across experiments?

%% HISTOGRAM OF ALL PAIRWISE CORRELATIONS
lt_figure; hold on;
title('distribution of pairwise correlations across all data');
xlabel('correlation coefficient');
ylabel('count (binsize: 0.1)');

Xcenters=(-1:0.1:0.9)+0.05;

% Plot all
Y = AllExptsData.pairwise_corr;
[hist_dat, ~] = hist(Y, Xcenters);

lt_plot_bar(Xcenters-0.01, hist_dat, {'Color', 'k', 'BarWidth', 0.4});

% Plot similar
inds=AllExptsData.pair_is_similar==1;
Y=AllExptsData.pairwise_corr(inds);
[hist_dat, ~] = hist(Y, Xcenters);

lt_plot_bar(Xcenters, hist_dat, {'Color', 'b', 'BarWidth', 0.4});

% Plot diff
inds=AllExptsData.pair_is_similar==0;
Y=AllExptsData.pairwise_corr(inds);
[hist_dat, ~] = hist(Y, Xcenters);

lt_plot_bar(Xcenters+0.01, hist_dat, {'Color', 'r', 'BarWidth', 0.4});

line([0 0], ylim, 'Color','k');


% ======= 1 pval vs. corr (all, similar, diff)


% ======== 2)

%% acoustic dist vs. song by song corr (ALL PAIRS)
hsplot=[];
lt_figure; hold on;

% similar
hsplot(1)=lt_subplot(3,3,1); hold on;
inds=AllExptsData.pair_is_similar==1;
X=AllExptsData.pairwise_corr(inds);
Y=AllExptsData.AcousticDist_all(inds);

lt_regress(Y, X, 1);
title('similar (all pairs)')
xlabel('corr (song by song)');
ylabel('acoustic euclidian dist');

% different
hsplot(2)=lt_subplot(3,3,2); hold on;
inds=AllExptsData.pair_is_similar==0;
X=AllExptsData.pairwise_corr(inds);
Y=AllExptsData.AcousticDist_all(inds);

lt_regress(Y, X, 1);
title('different (all pairs)');
xlabel('corr (song by song)');
ylabel('acoustic euclidian dist');

% all
hsplot(3)=lt_subplot(3,3,3); hold on;
inds=AllExptsData.pair_is_similar==1; % first simlilar
X=AllExptsData.pairwise_corr(inds);
Y=AllExptsData.AcousticDist_all(inds);

lt_plot(X, Y, {'Color', 'b'});

inds=AllExptsData.pair_is_similar==0; % first simlilar
X=AllExptsData.pairwise_corr(inds);
Y=AllExptsData.AcousticDist_all(inds);

lt_plot(X, Y, {'Color', 'r'});
title('all (all pairs)');
xlabel('corr (song by song)');
ylabel('acoustic euclidian dist');

% 3b) ---------- Use only within motif (motif by motif, song subtr)
% similar
hsplot(4)=lt_subplot(3,3,4); hold on;
inds=AllExptsData.pair_is_similar==1;
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_regress(Y, X, 1);
title('similar (same motif)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% similar
hsplot(5)=lt_subplot(3,3,5); hold on;
inds=AllExptsData.pair_is_similar==0;
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_regress(Y, X, 1);
title('different (same motif)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% all
hsplot(6)=lt_subplot(3,3,6); hold on;
title('all');
inds=AllExptsData.pair_is_similar==1; % first similar
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_plot(X, Y, {'Color', 'b'});

inds=AllExptsData.pair_is_similar==0; % second diff
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_plot(X, Y, {'Color', 'r'});
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% 3c - use absolute correlation instead of correlation.
hsplot(7)=lt_subplot(3,3,7); hold on;
inds=AllExptsData.pair_is_similar==1;
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);
X=abs(X);

lt_regress(Y, X, 1);
title('all (same motif, abs corr)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% similar
hsplot(8)=lt_subplot(3,3,8); hold on;
inds=AllExptsData.pair_is_similar==0;
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);
X=abs(X);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_regress(Y, X, 1);
title('different (same motif, abs corr)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% all
hsplot(9)=lt_subplot(3,3,9); hold on;
title('all');
inds=AllExptsData.pair_is_similar==1; % first similar
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);
X=abs(X);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_plot(X, Y, {'Color', 'b'});

inds=AllExptsData.pair_is_similar==0; % second diff
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);
X=abs(X);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_plot(X, Y, {'Color', 'r'});
title('all (same motif, abs corr)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% -------- FIGURE WIDE
% linkaxes(hsplot, 'xy');
lt_subtitle('Acoustic Dist vs. Corr (ignoring p-value)');


%% ================ ACOUSTIC DIST VS. CORR (only p<0.05)
hsplot=[];
lt_figure; hold on;

% ====== 1) first filter AllExptsData (only p<0.05)
% save old version to revert back to end of this code
AllExptsData_backup=AllExptsData;

indstokeep=AllExptsData.pairwise_corr_pval<0.05;

fieldsindstruct=fieldnames(AllExptsData);
for i=1:length(fieldsindstruct);
    fname=fieldsindstruct{i};
    AllExptsData.(fname)=AllExptsData.(fname)(indstokeep);
end



% ====== 2) Plot
% similar
hsplot(1)=lt_subplot(3,3,1); hold on;
inds=AllExptsData.pair_is_similar==1;
X=AllExptsData.pairwise_corr(inds);
Y=AllExptsData.AcousticDist_all(inds);

lt_regress(Y, X, 1);
title('similar (all pairs)')
xlabel('corr (song by song)');
ylabel('acoustic euclidian dist');

% different
hsplot(2)=lt_subplot(3,3,2); hold on;
inds=AllExptsData.pair_is_similar==0;
X=AllExptsData.pairwise_corr(inds);
Y=AllExptsData.AcousticDist_all(inds);

lt_regress(Y, X, 1);
title('different (all pairs)');
xlabel('corr (song by song)');
ylabel('acoustic euclidian dist');

% all
hsplot(3)=lt_subplot(3,3,3); hold on;
inds=AllExptsData.pair_is_similar==1; % first simlilar
X=AllExptsData.pairwise_corr(inds);
Y=AllExptsData.AcousticDist_all(inds);

lt_plot(X, Y, {'Color', 'b'});

inds=AllExptsData.pair_is_similar==0; % first simlilar
X=AllExptsData.pairwise_corr(inds);
Y=AllExptsData.AcousticDist_all(inds);

lt_plot(X, Y, {'Color', 'r'});
title('all (all pairs)');
xlabel('corr (song by song)');
ylabel('acoustic euclidian dist');

% 3b) ---------- Use only within motif (motif by motif, song subtr)
% similar
hsplot(4)=lt_subplot(3,3,4); hold on;
inds=AllExptsData.pair_is_similar==1;
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_regress(Y, X, 1);
title('similar (same motif)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% similar
hsplot(5)=lt_subplot(3,3,5); hold on;
inds=AllExptsData.pair_is_similar==0;
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_regress(Y, X, 1);
title('different (same motif)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% all
hsplot(6)=lt_subplot(3,3,6); hold on;
title('all');
inds=AllExptsData.pair_is_similar==1; % first similar
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_plot(X, Y, {'Color', 'b'});

inds=AllExptsData.pair_is_similar==0; % second diff
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_plot(X, Y, {'Color', 'r'});
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% 3c - use absolute correlation instead of correlation.
hsplot(7)=lt_subplot(3,3,7); hold on;
inds=AllExptsData.pair_is_similar==1;
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);
X=abs(X);

lt_regress(Y, X, 1);
title('all (same motif, abs corr)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% similar
hsplot(8)=lt_subplot(3,3,8); hold on;
inds=AllExptsData.pair_is_similar==0;
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);
X=abs(X);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_regress(Y, X, 1);
title('different (same motif, abs corr)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% all
hsplot(9)=lt_subplot(3,3,9); hold on;
title('all');
inds=AllExptsData.pair_is_similar==1; % first similar
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);
X=abs(X);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_plot(X, Y, {'Color', 'b'});

inds=AllExptsData.pair_is_similar==0; % second diff
X=AllExptsData.CorrCoeff_motif_SubtrSong_all(inds);
Y=AllExptsData.AcousticDist_all(inds);
X=abs(X);

inds2=~isnan(X);
X=X(inds2);
Y=Y(inds2);

lt_plot(X, Y, {'Color', 'r'});
title('all (same motif, abs corr)')
xlabel('corr (motif, song subtr)');
ylabel('acoustic euclidian dist');

% -------- FIGURE WIDE
% linkaxes(hsplot, 'xy');
lt_subtitle('Acoustic Dist vs. Corr (only p<0.05)');


% ====== 3) Revert back to old version
AllExptsData=AllExptsData_backup;


%% ============= 4) CORR VS. DISTANCE IN MOTIF


%% LEARNING vs. CORERLATION (1 for each expt, and 1 summary) (RELATIVE TO TARGET)

AllExptsData=struct;
AllExptsData.CorrCoeff_all=[];
AllExptsData.CorrCoeff_MotifMinusSong_all=[];
AllExptsData.Pval_all= [];
AllExptsData.Similar_status_all = [];
AllExptsData.AcoustDist_all=[];
AllExptsData.Pre_syl_sim_all=[];
AllExptsData.SylDistFromTarg=[];
AllExptsData.LearningRelTarg=[];

for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        SylFields_Unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        learning_targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        
        CorrCoeff_all=[];
        Pval_all=[];
        Similar_status_all=[];
        AcoustDist_all=[];
        Pre_syl_sim_all=[];
        SylDistFromTarg=[];
        LearningRelTarg=[];
        CorrCoeff_MotifMinusSong_all=[];
        
        for j=1:length(SylFields_Unique);
            syl=SylFields_Unique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                % skip is this is the target
                continue;
            end
            
            % ======= EXTRACT DATA
            % ---- Correlations
            corr=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(targsyl);
            p_val=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.p_val_vs.(targsyl);
            issimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            
            % --- corr using motifs
            try
                corr_motif_minussong=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(targsyl);
            catch err
                corr_motif_minussong=nan;
            end
            
            
            % ---- Other things.
            acoust_dist_Z=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
            pre_syl_similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            distance_from_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ;
            
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learning_rel_targ=learning/learning_targsyl;

            
            
            % ===================  OUTPUT
            if ~isempty(corr);
                CorrCoeff_all=[CorrCoeff_all corr];
                CorrCoeff_MotifMinusSong_all=[CorrCoeff_MotifMinusSong_all corr_motif_minussong];
                Pval_all=[Pval_all p_val];
                Similar_status_all=[Similar_status_all issimilar];
                AcoustDist_all=[AcoustDist_all acoust_dist_Z];
                Pre_syl_sim_all=[Pre_syl_sim_all pre_syl_similar];
                SylDistFromTarg=[SylDistFromTarg distance_from_targ];
                LearningRelTarg=[LearningRelTarg learning_rel_targ];
                
                
                
                % ===== OUTPUT FOR ALL DATA ================
                AllExptsData.CorrCoeff_all=[AllExptsData.CorrCoeff_all corr];
                AllExptsData.CorrCoeff_MotifMinusSong_all=[AllExptsData.CorrCoeff_MotifMinusSong_all corr_motif_minussong];
                AllExptsData.Pval_all= [AllExptsData.Pval_all p_val];
                AllExptsData.Similar_status_all = [AllExptsData.Similar_status_all issimilar];
                AllExptsData.AcoustDist_all=[AllExptsData.AcoustDist_all acoust_dist_Z];
                AllExptsData.Pre_syl_sim_all=[AllExptsData.Pre_syl_sim_all pre_syl_similar];
                AllExptsData.SylDistFromTarg=[AllExptsData.SylDistFromTarg distance_from_targ];
                AllExptsData.LearningRelTarg=[AllExptsData.LearningRelTarg learning_rel_targ];
                
                % =========================================
            end
            
        end
        
        % ====== PLOT (for this bird/epxt)
        lt_figure; hold on;
        title([birdname ' ' exptname]);
        
        % 1) ------ Corr coeff vs. learning.
        lt_subplot(1,3,1); hold on; title('learning vs. corr coeff');
        xlabel('corr coeff');
        ylabel('learning rel. target');
        
        % similar
        inds=Similar_status_all==1;
        X = CorrCoeff_all(inds);
        Y = LearningRelTarg(inds);
        
        if ~isempty(X);
            lt_plot(X, Y, {'Color' ,'b'});
        end
        
        % diff
        inds=Similar_status_all==0;
        X = CorrCoeff_all(inds);
        Y = LearningRelTarg(inds);
        
        if ~isempty(X);
            lt_plot(X, Y, {'Color' ,'r'});
        end
        
        line([0 0], ylim);
        line(xlim, [0 0]);
        
        
        % 2) ------- Corr coeff vs. acoustic distance
        lt_subplot(1,3,2); hold on; title('acoustic distance vs. corr coeff');
        xlabel('corr coeff');
        ylabel('acoustic distance');
        
        % similar
        inds=Similar_status_all==1;
        X = CorrCoeff_all(inds);
        Y = AcoustDist_all(inds);
        
        if ~isempty(X);
            lt_plot(X, Y, {'Color' ,'b'});
        end
        
        % diff
        inds=Similar_status_all==0;
        X = CorrCoeff_all(inds);
        Y = AcoustDist_all(inds);
        
        if ~isempty(X);
            lt_plot(X, Y, {'Color' ,'r'});
        end
        
        line([0 0], ylim);
        line(xlim, [0 0]);
        
        lt_subtitle([birdname ' ' exptname]);
        
        
        
    end
end




%% ====================== PLOT FOR ALL DATA (all expts, all birds)
only_significant_corr_SummaryPlots=0; % then only corrs with p<0.05;
alpha=1;

lt_figure; hold on;
hsplot=[];

% ======= ALL DATA
% 1) CORR (song) vs. learning.
hsplot(1)=lt_subplot(3,3,1); hold on; title('all');
xlabel('corr coeff');
ylabel('learning rel. target');

X = AllExptsData.CorrCoeff_all;
Y = AllExptsData.LearningRelTarg;
pvals=AllExptsData.Pval_all;
similar=AllExptsData.Similar_status_all;


% only significnat correlatiosn?
if only_significant_corr_SummaryPlots==1;
    inds=pvals<alpha;
    X=X(inds);
    Y=Y(inds);
    similar=similar(inds);
end
    
% SUBPLOT 1) all
lt_regress(Y,X,1);

% similar
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end

% diff
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end

line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 2) similar
% similar
lt_subplot(3,3,2); hold on;
title('similar syls');
xlabel('corr (song)');
ylabel('learning');
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 3) different
lt_subplot(3,3,3); hold on;
title('different syls');
xlabel('corr (song)');
ylabel('learning');
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% ===== USING MOTIF CORRELATIONS
lt_subplot(3,3,4); hold on; title('all');
xlabel('corr coeff (motif-song)');
ylabel('learning rel. target');

X = AllExptsData.CorrCoeff_MotifMinusSong_all;
Y = AllExptsData.LearningRelTarg;
pvals=AllExptsData.Pval_all;
similar=AllExptsData.Similar_status_all;


% only significnat correlatiosn?
if only_significant_corr_SummaryPlots==1;
    inds=pvals<alpha;
    X=X(inds);
    Y=Y(inds);
    similar=similar(inds);
end
    
% SUBPLOT 1) all
lt_regress(Y,X,1);

% similar
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end

% diff
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end

line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 2) similar
% similar
lt_subplot(3,3,5); hold on;
title('similar syls');
xlabel('corr (motif-song)');
ylabel('learning');
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 3) different
lt_subplot(3,3,6); hold on;
title('different syls');
xlabel('corr (motif-song)');
ylabel('learning');
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% ==== ONLY SYLS ON DIFFERENT MOTIF (SONG CORR)
lt_subplot(3,3,7); hold on; title('all syls(diff motif)');
xlabel('corr coeff');
ylabel('learning rel. target');

X = AllExptsData.CorrCoeff_all;
Y = AllExptsData.LearningRelTarg;
pvals=AllExptsData.Pval_all;
motif=AllExptsData.SylDistFromTarg;
similar=AllExptsData.Similar_status_all;

% only those on different motif
inds=isnan(motif);
X=X(inds);
Y=Y(inds);
similar=similar(inds);
pvals=pvals(inds);

% only significnat correlatiosn?
if only_significant_corr_SummaryPlots==1;
    inds=pvals<alpha;
    X=X(inds);
    Y=Y(inds);
    similar=similar(inds);
end
    
% SUBPLOT 1) all
lt_regress(Y,X,1);

% similar
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end

% diff
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end

line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 2) similar
% similar
lt_subplot(3,3,8); hold on;
title('similar syls (diffmotif)');
xlabel('corr (song)');
ylabel('learning');
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 3) different
lt_subplot(3,3,9); hold on;
title('different syls (diffmotif)');
xlabel('corr (song)');
ylabel('learning');
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% === figurewide
if only_significant_corr_SummaryPlots==1;
    lt_subtitle(['learning vs. corr, only corrs with p< ' num2str(alpha) ]);
else
    lt_subtitle('learning vs. corr, dont care about p-value');
end

%% ====================== PLOT FOR ALL DATA (all expts, all birds) [ONLY SIGNIFICANT CORRS]
only_significant_corr_SummaryPlots=1; % then only corrs with p<0.05;
alpha=0.05;

lt_figure; hold on;
hsplot=[];

% ======= ALL DATA
% 1) CORR (song) vs. learning.
hsplot(1)=lt_subplot(3,3,1); hold on; title('all');
xlabel('corr coeff');
ylabel('learning rel. target');

X = AllExptsData.CorrCoeff_all;
Y = AllExptsData.LearningRelTarg;
pvals=AllExptsData.Pval_all;
similar=AllExptsData.Similar_status_all;


% only significnat correlatiosn?
if only_significant_corr_SummaryPlots==1;
    inds=pvals<alpha;
    X=X(inds);
    Y=Y(inds);
    similar=similar(inds);
end
    
% SUBPLOT 1) all
lt_regress(Y,X,1);

% similar
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end

% diff
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end

line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 2) similar
% similar
lt_subplot(3,3,2); hold on;
title('similar syls');
xlabel('corr (song)');
ylabel('learning');
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 3) different
lt_subplot(3,3,3); hold on;
title('different syls');
xlabel('corr (song)');
ylabel('learning');
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% ===== USING MOTIF CORRELATIONS
lt_subplot(3,3,4); hold on; title('all');
xlabel('corr coeff (motif-song)');
ylabel('learning rel. target');

X = AllExptsData.CorrCoeff_MotifMinusSong_all;
Y = AllExptsData.LearningRelTarg;
pvals=AllExptsData.Pval_all;
similar=AllExptsData.Similar_status_all;


% only significnat correlatiosn?
if only_significant_corr_SummaryPlots==1;
    inds=pvals<alpha;
    X=X(inds);
    Y=Y(inds);
    similar=similar(inds);
end
    
% SUBPLOT 1) all
lt_regress(Y,X,1);

% similar
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end

% diff
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end

line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 2) similar
% similar
lt_subplot(3,3,5); hold on;
title('similar syls');
xlabel('corr (motif-song)');
ylabel('learning');
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 3) different
lt_subplot(3,3,6); hold on;
title('different syls');
xlabel('corr (motif-song)');
ylabel('learning');
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% ==== ONLY SYLS ON DIFFERENT MOTIF (SONG CORR)
lt_subplot(3,3,7); hold on; title('all syls(diff motif)');
xlabel('corr coeff');
ylabel('learning rel. target');

X = AllExptsData.CorrCoeff_all;
Y = AllExptsData.LearningRelTarg;
pvals=AllExptsData.Pval_all;
motif=AllExptsData.SylDistFromTarg;
similar=AllExptsData.Similar_status_all;

% only those on different motif
inds=isnan(motif);
X=X(inds);
Y=Y(inds);
similar=similar(inds);
pvals=pvals(inds);

% only significnat correlatiosn?
if only_significant_corr_SummaryPlots==1;
    inds=pvals<alpha;
    X=X(inds);
    Y=Y(inds);
    similar=similar(inds);
end
    
% SUBPLOT 1) all
lt_regress(Y,X,1);

% similar
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end

% diff
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end

line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 2) similar
% similar
lt_subplot(3,3,8); hold on;
title('similar syls (diffmotif)');
xlabel('corr (song)');
ylabel('learning');
inds=similar==1;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'b'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% SUBPLOT 3) different
lt_subplot(3,3,9); hold on;
title('different syls (diffmotif)');
xlabel('corr (song)');
ylabel('learning');
inds=similar==0;
Xtmp = X(inds);
Ytmp = Y(inds);

if ~isempty(Xtmp);
    lt_regress(Ytmp,Xtmp,1);
    lt_plot(Xtmp, Ytmp, {'Color' ,'r'});
end
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');

% === figurewide
if only_significant_corr_SummaryPlots==1;
    lt_subtitle(['learning vs. corr, only corrs with p< ' num2str(alpha) ]);
else
    lt_subtitle('learning vs. corr, dont care about p-value');
end

%% (VS TARGET) OVER ALL EXPERIMENTS - PLOT LEARNING VS. COM DISTANCE FROM TARGET SYL
write_syl=0; % if 1, then writes syl name (with target, bird,e tc)

% USING CONSOLID START FOR LEARNING
xlabel('eucl acoustic distance');
ylabel('learning (rel to target)');

Learning_all=[];
Acoustic_all=[];
Similar_all=[];
BirdNums_all=[];
Targets_all={};
Syls_all={};
ExptNum_all=[];

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        % ==== targsyl stats
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
        
        learning_targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        
        % ==== PLOT ALL SYLS
        syllist=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed);
        
        for j=1:length(syllist);
            syl=syllist{j};
            
            % skip the target
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue;
            end
            
            
            % ==== collect data
            % - Learning
            learning_notreltarg=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learning=learning_notreltarg/learning_targsyl;

            Learning_all=[Learning_all learning];
            
            % - Distance
            eucl_dist=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
            Acoustic_all=[Acoustic_all eucl_dist];
            
            % is similar?
            is_similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            Similar_all=[Similar_all is_similar];
            
            
            % Things for writing
            BirdNums_all=[BirdNums_all i];
            Targets_all=[Targets_all {targsyl}];
            Syls_all=[Syls_all {syl}];
            ExptNum_all=[ExptNum_all, ii];
            
            
        end
    end
end
lt_plot_zeroline;

% ==== PLOT 
lt_figure; hold on;

% -- ALL
lt_subplot(1,3,1); hold on; grid on;
title('All');
xlabel('acoust dist');
ylabel('learning metric');

X=Acoustic_all;
Y=Learning_all;

lt_regress(Y,X,1);

% plot names 
if write_syl==1;
   for j=1:length(BirdNums_all);
       string=[num2str(BirdNums_all(j)) ' ' num2str(ExptNum_all(j)) ' ' Syls_all{j} ' ' Targets_all{j}];
       
       text(X(j)+0.1, Y(j), string);
   end
    
    
end
% -- SIMILAR
lt_subplot(1,3,2); hold on; grid on;
title('Similar');
xlabel('acoust dist');
ylabel('learning');

inds=Similar_all==1;
X=Acoustic_all(inds);
Y=Learning_all(inds);

lt_regress(Y,X,1);
lt_plot(X,Y, {'Color', 'b'});

subplot(1,3,1); 
lt_plot(X,Y, {'Color', 'b'});

% -- DIFFERENT
lt_subplot(1,3,3); hold on; grid on;
title('Different');
xlabel('acoust dist');
ylabel('learning');

inds=Similar_all==0;
X=Acoustic_all(inds);
Y=Learning_all(inds);

lt_regress(Y,X,1);
lt_plot(X,Y, {'Color', 'r'});
subplot(1,3,1); 
lt_plot(X,Y, {'Color', 'r'});

% ----
lt_subtitle('Learning vs. acoustic distance');


%% PLOT DISTRIBUTION OF ACOUSTIC SIMILARITIES

lt_figure; hold on;

BinCenters=linspace(min(Acoustic_all)-0.2, max(Acoustic_all)+0.2, 20);

lt_subplot(1,2,1); hold on;
title('similar and diff');
% plot similar
inds=Similar_all==1;
Y=Acoustic_all(inds);

Yhist=hist(Y, BinCenters);

lt_plot_bar(BinCenters, Yhist, {'Color','b'});

% plot different
inds=Similar_all==0;
Y=Acoustic_all(inds);

Yhist=hist(Y, BinCenters);

lt_plot_bar(BinCenters, Yhist, {'Color','r'});

% plot all
lt_subplot(1,2,2); hold on; title('all');

Y=Acoustic_all;

Yhist=hist(Y, BinCenters);

lt_plot_bar(BinCenters, Yhist, {'Color','k'});

lt_subtitle('All acoustic distances (rel target)');


%% PLOT LEARNING REL TO DISTANCE IN MOTIF

write_corr=0;
write_syl=0;

lt_figure; hold on;
title('All syls, all expts (consol start mean FF)');
xlabel('acoustic distance');
ylabel('learning (rel to target)');

Learning_AllExpts=[];
MotifDist_AllExpts=[];
Similar_AllExpts=[];
PreSylSimilar_AllExpts=[];

lt_subplot(2,4,1); hold on; grid on
title('not in target motif')
xlim([-1 1]);
ylim([-1 1]);
ylabel('learning (rel to target)');

lt_subplot(2,4,2); hold on; grid on
title('in target motif')
xlim([-1 1]);
ylim([-1 1]);
xlabel('motif distance');
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
xlabel('motif distance');
ylabel('learning (rel to target)');

lt_subplot(2,4,7:8); hold on; grid on
title('in target motif');
ylim([-0.5 0.5]);
xlabel('motif distance');
ylabel('learning (rel to target)');

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        learning_targ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
        
        % ==== PLOT ALL SYLS
        syllist=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed);
        
        for j=1:length(syllist);
            syl=syllist{j};
            
            % skip if this is a target
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue;
            end
            
            % ==== collect data
            learning=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            learning_rel_targ=learning/learning_targ;

            if isnan(learning_rel_targ);
                keyboard;
                dascadcv
            end
            
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
learning_mean=mean(learning);
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

Ymean=mean(Y);
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

end

%% === STATISTICAL TESTS ( ON LAERNING VS. MOTIF STUFF) (1. in vs. out of motif) (2. POSITION IN MOTIF)
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

% group_anova{2} = Similar_AllExpts(inds); % similar or different
group_anova{1} = MotifDist_AllExpts(inds); % position on motif (or on different motif)
Y = Learning_AllExpts(inds);

GroupNames={'similarity', 'motifpos'};
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



