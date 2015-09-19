function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMAN(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 7/21/15 - Plots LMAN learning
close all;

%% SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA

% copy strcuture, save backup.
if ~exist('SeqDepPitch_AcrossBirds_ORIGINAL', 'var'); % do not save a new backup if it already exists.
    SeqDepPitch_AcrossBirds_ORIGINAL=SeqDepPitch_AcrossBirds;
end

filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);



%% CORRELATIONS - COMPARE MUSC VS. PBS (BAR PLOTS)

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
            
            
            % ========== PBS DATA
            Y_corr=[];
            X=[];
            
            % ---- GET CORR AGAINST OTHER SYLS (SONG BY SONG)
            Syllist_all=[];
            motif_divider=[];
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                X=[X jj];
                Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(othersyl)];
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(othersyl));
                    Y_corr=[Y_corr nan];
                end
                
                Syllist_all=[Syllist_all {SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).single_syl}];
            end
            
            % PLOT
            lt_plot_bar(X-0.3, Y_corr, {'Color', 'w', 'BarWidth', 0.25});
            
            
            % ========== MUSC DATA
            Y_corr=[];
            X=[];
            
            % ---- GET CORR AGAINST OTHER SYLS (SONG BY SONG)
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                X=[X jj];
                Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(othersyl)];
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(othersyl));
                    Y_corr=[Y_corr nan];
                end
            end
            
            % PLOT
            lt_plot_bar(X-0.1, Y_corr,{'Color', 'k', 'BarWidth', 0.25});
            
            ylim([-0.4, 1]);
            set(gca, 'XTick', 1:length(Syllist_all))
            set(gca, 'XTickLabel', Syllist_all);
            
            
            % ============== GET CORR WITHIN MOTIF (PBS)
            X=[];
            Y_corr=[];
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                
                % make sure this syl is in the same motif
                
                X=[X jj];
                try
                    Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(othersyl)];
                catch err
                    Y_corr=[Y_corr nan];
                end
                
            end
            
            lt_plot_bar(X+0.1, Y_corr, {'Color', 'y', 'BarWidth', 0.25});
            
            
            % ============== GET CORR WITHIN MOTIF (MUSC)
            X=[];
            Y_corr=[];
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                
                % make sure this syl is in the same motif
                
                X=[X jj];
                try
                    Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(othersyl)];
                catch err
                    Y_corr=[Y_corr nan];
                end
                
            end
            
            lt_plot_bar(X+0.3, Y_corr, {'Color', 'b', 'BarWidth', 0.25});
            
        end
        
        lt_subtitle(['Correlations (LtoR: song(PBS, MUSC), motif(PBS, MUSC)); Bird: ' birdname ', expt ' exptname]);
    end
end


%% CORRELATION P VALUES - COMPARE MUSC VS. PBS (BAR PLOTS)

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
            
            
            % ========== PBS DATA
            Y_corr=[];
            X=[];
            
            % ---- GET CORR AGAINST OTHER SYLS (SONG BY SONG)
            Syllist_all=[];
            motif_divider=[];
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                X=[X jj];
                Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.p_val_vs.(othersyl)];
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.p_val_vs.(othersyl));
                    Y_corr=[Y_corr nan];
                end
                
                Syllist_all=[Syllist_all {SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(othersyl).single_syl}];
            end
            
            % PLOT
            lt_plot_bar(X-0.3, log10(Y_corr), {'Color', 'w', 'BarWidth', 0.25});
            
            
            % ========== MUSC DATA
            Y_corr=[];
            X=[];
            
            % ---- GET CORR AGAINST OTHER SYLS (SONG BY SONG)
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                X=[X jj];
                Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.p_val_vs.(othersyl)];
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.song_by_song.p_val_vs.(othersyl));
                    Y_corr=[Y_corr nan];
                end
            end
            
            % PLOT
            lt_plot_bar(X-0.1, log10(Y_corr),{'Color', 'k', 'BarWidth', 0.25});
            
            %             ylim([-0.4, 1]);
            set(gca, 'XTick', 1:length(Syllist_all))
            set(gca, 'XTickLabel', Syllist_all);
            
            
            % ============== GET CORR WITHIN MOTIF (PBS)
            X=[];
            Y_corr=[];
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                
                % make sure this syl is in the same motif
                
                X=[X jj];
                try
                    Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.p_val_vs.(othersyl)];
                catch err
                    Y_corr=[Y_corr nan];
                end
                
            end
            
            lt_plot_bar(X+0.1, log10(Y_corr), {'Color', 'y', 'BarWidth', 0.25});
            
            
            % ============== GET CORR WITHIN MOTIF (MUSC)
            X=[];
            Y_corr=[];
            for jj=1:length(SylFields_Unique);
                othersyl=SylFields_Unique{jj};
                
                % make sure this syl is in the same motif
                
                X=[X jj];
                try
                    Y_corr=[Y_corr SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions_LMAN.(syl).CORRELATIONS.motif_by_motif.p_val_vs.(othersyl)];
                catch err
                    Y_corr=[Y_corr nan];
                end
                
            end
            
            lt_plot_bar(X+0.3, log10(Y_corr), {'Color', 'b', 'BarWidth', 0.25});
            
            
            % Line for p=0.05
            line(xlim, [log10(0.05) log10(0.05)])
        end
        
        lt_subtitle(['corr p-values (LtoR: song(PBS, MUSC), motif(PBS, MUSC)); Bird: ' birdname ', expt ' exptname]);
    end
end



%% COMPARE ALL PAIRWISE CORRELATIONS

        conditions={'' , '_LMAN'}; % pbs/musc

% AllExptsData.AcousticDist_all=[];

AllExptsData=struct;
AllExptsData.pairwise_corr=cell(length(conditions), 1);
AllExptsData.pairwise_corr_pval=cell(length(conditions), 1);
AllExptsData.pair_is_similar=cell(length(conditions), 1);
AllExptsData.CorrCoeff_motif_all=cell(length(conditions), 1);
AllExptsData.Pval_motif_all=cell(length(conditions), 1);
AllExptsData.CorrCoeff_motif_SubtrSong_all=cell(length(conditions), 1);
AllExptsData.Pval_motif_SubtrSong_all=cell(length(conditions), 1);

AllExptsData.bird_expt=cell(length(conditions), 1);

% Initiate figures that plot all birds, each expt one point
hfig_dots=lt_figure(); hold on;


for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        SylFields_Unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % to hold data for this experiment (cell to hold either pbs(1) or
        % musc(2)
        CorrCoeff_all={};
        Pval_all={};
        Similar_status_all={};
        %             AcousticDist_all={};
        CorrCoeff_motif_all={};
        Pval_motif_all={};
        CorrCoeff_motif_SubtrSong_all={};
        Pval_motif_SubtrSong_all={};
        
        % --- PBS or MUSC
        for iii=1:length(conditions)
            sylIDfield=['Syl_ID_Dimensions' conditions{iii}]; % pbs or musc
            
            if isempty(AllExptsData.pairwise_corr{iii}); % only do this once
                AllExptsData.pairwise_corr{iii}=[];
                AllExptsData.pairwise_corr_pval{iii}=[];
                AllExptsData.pair_is_similar{iii}=[];
                AllExptsData.CorrCoeff_motif_all{iii}=[];
                AllExptsData.Pval_motif_all{iii}=[];
                AllExptsData.CorrCoeff_motif_SubtrSong_all{iii}=[];
                AllExptsData.Pval_motif_SubtrSong_all{iii}=[];
                AllExptsData.bird_expt{iii}=[];
            end
            
            CorrCoeff_all{iii}=[];
            Pval_all{iii}=[];
            Similar_status_all{iii}=[];
            
            CorrCoeff_motif_all{iii}=[];
            Pval_motif_all{iii}=[];
            
            CorrCoeff_motif_SubtrSong_all{iii}=[];
            Pval_motif_SubtrSong_all{iii}=[];
            
            
            % DATA FOR THIS EXPERIMENT -------------------------------------------------------------------
            % SYLLABLES
            for j=1:length(SylFields_Unique);
                syl=SylFields_Unique{j};
                
                % =========  GET CORR AGAINST ALL SYLLABLES THAT COME AFTER IT (no double
                % counting, and no auto-corr).
                if j==length(SylFields_Unique); % last syl - stop.
                    continue;
                end
                
                % ---- OTHER SYLLABLES
                for jj=j+1:length(SylFields_Unique);
                    othersyl=SylFields_Unique{jj};
                    
                    % === EXTRACT DATA
                    % -- Corr (song-by_song)
                    corr=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(othersyl);
                    p_val=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(syl).CORRELATIONS.song_by_song.p_val_vs.(othersyl);
                    
                    % --- Corr (motif_by_motif)
                    if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs, othersyl); % might not be in same motif
                        corr_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(othersyl);
                        pval_motif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(syl).CORRELATIONS.motif_by_motif.p_val_vs.(othersyl);
                        
                        corr_motif_MinusSongMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(othersyl);
                        pval_motif_MinusSongMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.p_val_vs.(othersyl);
                    else
                        corr_motif=nan;
                        pval_motif=nan;
                        
                        corr_motif_MinusSongMean=nan;
                        pval_motif_MinusSongMean=nan;
                    end
                    
                    % ---- Similar vs different
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(syl).single_syl == SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.(sylIDfield).(othersyl).single_syl
                        issimilar=1;
                    else
                        issimilar=0;
                    end
                    
                    
                    % ====== OUTPUT
                    if ~isempty(corr) && ~isnan(corr);
                        CorrCoeff_all{iii}=[CorrCoeff_all{iii} corr];
                        Pval_all{iii}=[Pval_all{iii} p_val];
                        Similar_status_all{iii}=[Similar_status_all{iii} issimilar];
                        %                     AcousticDist_all=[AcousticDist_all acoustic_dist];
                        
                        CorrCoeff_motif_all{iii}=[CorrCoeff_motif_all{iii} corr_motif];
                        Pval_motif_all{iii}=[Pval_motif_all{iii} pval_motif];
                        
                        CorrCoeff_motif_SubtrSong_all{iii}=[CorrCoeff_motif_SubtrSong_all{iii} corr_motif_MinusSongMean];
                        Pval_motif_SubtrSong_all{iii}=[Pval_motif_SubtrSong_all{iii} pval_motif_MinusSongMean];
                        

                    end
                end
            end
            
        end
        
        % ==================== PUT INTO ALL EXPTS STRUCTURE (one ind for
        % MUSC PBS each) (across birds/expts)
        for iii=1:length(conditions)
        AllExptsData.pairwise_corr{iii}=[AllExptsData.pairwise_corr{iii} CorrCoeff_all{iii}];
        AllExptsData.pairwise_corr_pval{iii}=[AllExptsData.pairwise_corr_pval{iii} Pval_all{iii}];
        
        AllExptsData.pair_is_similar{iii}=[AllExptsData.pair_is_similar{iii} Similar_status_all{iii}];
        
        AllExptsData.CorrCoeff_motif_all{iii}=[AllExptsData.CorrCoeff_motif_all{iii} CorrCoeff_motif_all{iii}];
        AllExptsData.Pval_motif_all{iii}=[AllExptsData.Pval_motif_all{iii} Pval_motif_all{iii}];
        
        AllExptsData.CorrCoeff_motif_SubtrSong_all{iii}=[AllExptsData.CorrCoeff_motif_SubtrSong_all{iii} CorrCoeff_motif_SubtrSong_all{iii}];
        AllExptsData.Pval_motif_SubtrSong_all{iii}=[AllExptsData.Pval_motif_SubtrSong_all{iii} Pval_motif_SubtrSong_all{iii}];
        AllExptsData.bird_expt{iii}=[AllExptsData.bird_expt{iii}; repmat([i ii], length(Pval_motif_SubtrSong_all{iii}),1)];
        end
        
        % =============== PLOT 1) SCATTER, all vals (song) for this bird
        lt_figure; hold on;
        Xcorr=CorrCoeff_all{1};
        Ycorr=CorrCoeff_all{2};
        
        % ------------------ SCATTER
        % --- all
        lt_subplot(2,3,1); hold on; box on; grid on;
        title('all pairs');
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        
        lt_plot(X, Y);
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        
        % --- similar
        lt_subplot(2,3,2); hold on; box on; grid on;
        title('similar pairs');
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        
        similar=Similar_status_all{iii};
        inds=similar==1;
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y, {'Color','b'});
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        
        
        % --- different
        lt_subplot(2,3,3); hold on; box on; grid on;
        title('diff pairs');
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        
        similar=Similar_status_all{iii};
        inds=similar==0;
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y, {'Color','r'});
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        
        
        % ------------------ PLOT Only those passing p=alpha threshold
        % --- all
        alpha=0.005; % to filter pairs
        lt_subplot(2,3,4); hold on; box on; grid on;
        title(['all pairs (p< ' num2str(alpha) ')']);
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        inds=(Pval_all{1}<alpha & Pval_all{2}<alpha);
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y);
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        % --- similar
        lt_subplot(2,3,5); hold on; box on; grid on;
        title(['similar pairs (p< ' num2str(alpha) ')']);
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        similar=Similar_status_all{iii};
        
        inds=(Pval_all{1}<alpha & Pval_all{2}<alpha); % filter p value
        X=X(inds);
        Y=Y(inds);
        similar=similar(inds);
        
        inds=similar==1;
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y, {'Color','b'});
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        
        % --- different
        lt_subplot(2,3,6); hold on; box on; grid on;
        title(['diff pairs (p< ' num2str(alpha) ')']);
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        similar=Similar_status_all{iii};
        
        inds=(Pval_all{1}<alpha & Pval_all{2}<alpha); % filter p value
        X=X(inds);
        Y=Y(inds);
        similar=similar(inds);
        
        
        inds=similar==0;
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y, {'Color','r'});
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        % --- FIGUREWIDE
        lt_subtitle(['Song corr, ' birdname '-' exptname]);
        
        
        % ================ PLOT 2) scatter, using motif, for this bird
        lt_figure; hold on;
        Xcorr=CorrCoeff_motif_all{1};
        Ycorr=CorrCoeff_motif_all{2};
        
        % ------------------ SCATTER
        % --- all
        lt_subplot(2,3,1); hold on; box on; grid on;
        title('all pairs');
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        
        lt_plot(X, Y);
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        % --- similar
        lt_subplot(2,3,2); hold on; box on; grid on;
        title('similar pairs');
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        
        similar=Similar_status_all{iii};
        inds=similar==1;
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y, {'Color','b'});
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        
        % --- different
        lt_subplot(2,3,3); hold on; box on; grid on;
        title('diff pairs');
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        
        similar=Similar_status_all{iii};
        inds=similar==0;
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y, {'Color','r'});
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        
        % --------------------- PLOT Only those passing p=alpha threshold
        % --- all
        alpha=0.005; % to filter pairs
        lt_subplot(2,3,4); hold on; box on; grid on;
        title(['all pairs (p< ' num2str(alpha) ')']);
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        inds=(Pval_all{1}<alpha & Pval_all{2}<alpha);
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y);
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        % --- similar
        lt_subplot(2,3,5); hold on; box on; grid on;
        title(['similar pairs (p< ' num2str(alpha) ')']);
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        similar=Similar_status_all{iii};
        
        inds=(Pval_all{1}<alpha & Pval_all{2}<alpha); % filter p value
        X=X(inds);
        Y=Y(inds);
        similar=similar(inds);
        
        inds=similar==1;
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y, {'Color','b'});
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        
        % --- different
        lt_subplot(2,3,6); hold on; box on; grid on;
        title(['diff pairs (p< ' num2str(alpha) ')']);
        xlabel('corr (song) (PBS)');
        ylabel('corr (song) (MUSC)');
        
        X=Xcorr;
        Y=Ycorr;
        similar=Similar_status_all{iii};
        
        inds=(Pval_all{1}<alpha & Pval_all{2}<alpha); % filter p value
        X=X(inds);
        Y=Y(inds);
        similar=similar(inds);
        
        
        inds=similar==0;
        X=X(inds);
        Y=Y(inds);
        
        lt_plot(X, Y, {'Color','r'});
        xlim([-0.5 1]);
        ylim([-0.5 1]);
        
        line([-0.5 1],[-0.5 1]);
        
        % p value
        try
        p=signrank(X,Y);
        lt_plot_pvalue(p);
        catch err
        end
        
        % --- FIGUREWIDE
        lt_subtitle(['Motif corr, ' birdname '-' exptname]);
        
        
        % ============== PLOT 2) One value for each bird
        figure(hfig_dots); hold on;
        Xcorr=CorrCoeff_all{1};
        Ycorr=CorrCoeff_all{2};
        
        % ------------------ SCATTER
        % --- all
        lt_subplot(2,3,1); hold on; box on; grid on;
        title('all pairs');
        ylabel('corr (song)');
        
        X=Xcorr;
        Y=Ycorr;
        
        Xmean=nanmean(X);
        Xsem=lt_sem(X);
        
        Ymean=nanmean(Y);
        Ysem=lt_sem(Y);
        
        errorbar([1 2], [Xmean Ymean], [Xsem Ysem], 'o-', 'MarkerSize', 8);
        
        % p value
        try
        p=signrank(X,Y);
%         lt_plot_pvalue(p);
        catch err
        end

        text(2.1, Ymean, [birdname(1:4) '-' exptname(end-8:end) ', p=' num2str(p)], 'FontSize', 12', 'FontWeight', 'bold');
        xlim([0 4]);
        ylim([-0.5 1]);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'PBS', 'MUSC'})
        
        lt_plot_zeroline;
        
        
%         % --- similar
%         lt_subplot(2,3,1); hold on; box on; grid on;
%         title('all pairs');
%         ylabel('corr (song)');
%         
%         X=Xcorr;
%         Y=Ycorr;
%         
%         Xmean=nanmean(X);
%         Xsem=lt_sem(X);
%         
%         Ymean=nanmean(Y);
%         Ysem=lt_sem(Y);
%         
%         errorbar([1 2], [Xmean Ymean], [Xsem Ysem], 'o-', 'MarkerSize', 8);
%         
%                 p=signrank(X,Y);
% %         lt_plot_pvalue(p);
% 
%         text(2.1, [birdname '-' exptname 'p=' num2str(p)], 'FontSize', 12', 'FontWeight', 'bold');
% 
%         
%         xlim([0 4]);
%         ylim([-0.5 1]);
%                 set(gca, 'XTick', [1 2], 'XTickLabel', {'PBS', 'MUSC'})
% 
%         lt_plot_zeroline;
        
        % ===== PLOT Only those passing p=alpha threshold
        % --- all
        alpha=0.005; % to filter pairs
        lt_subplot(2,3,4); hold on; box on; grid on;
        title(['all pairs (p< ' num2str(alpha) ')']);
        ylabel('corr (song)');
        
        X=Xcorr;
        Y=Ycorr;
        inds=(Pval_all{1}<alpha & Pval_all{2}<alpha);
        X=X(inds);
        Y=Y(inds);
        
        Xmean=nanmean(X);
        Xsem=lt_sem(X);
        
        Ymean=nanmean(Y);
        Ysem=lt_sem(Y);
        
        errorbar([1 2], [Xmean Ymean], [Xsem Ysem], 'o-', 'MarkerSize', 8);
        
        try
        p=signrank(X,Y);
%         lt_plot_pvalue(p);
        catch err
        end

        text(2.1, Ymean,  [birdname(1:4) '-' exptname(end-8:end) ', p=' num2str(p)], 'FontSize', 12', 'FontWeight', 'bold');

        xlim([0 4]);
        ylim([-0.5 1]);
        
        lt_plot_zeroline;
        

        % FIGUREWIDE
        lt_subtitle('PBS vs. MUSC - Each bird/expt one mean+sem');
    end
end


% ========= PLOT ACROSS BIRDS
% PLOT 1) One value (mean, sem) per expt


%% PLOT vs. motif position