function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_MULTIDIR(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 7/19/15 - Plots multidir learning experiments (starting multidir or adding target during experiment)

%% EXTRACT ONLY EXPERIMENTS WITH MULTIDIR
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

expts_to_remove=[];
for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    expts_to_remove{i}=[];
    
    for ii=1:NumExperiments
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart');
            expts_to_remove{i}=[expts_to_remove{i} ii];
        end
    end
end

% actually remove those experiments from structure
for i=1:length(expts_to_remove);
    if ~isempty(expts_to_remove{i});
        SeqDepPitch_AcrossBirds.birds{i}.experiment(expts_to_remove{i})=[];
        disp(['TO GET ONLY MULTIDIR EXPERIMENTS, REMOVED: ' num2str(i) '; expt: ' num2str(expts_to_remove{i})]);
    end
end


PARAMS.global.MULTIDIR.expts_to_remove=expts_to_remove;
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

% total num experiments
TotalNumExperiments=0;
for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:numexperiments;
        TotalNumExperiments=TotalNumExperiments+1;
    end
end


%% PLOT ALL EXPERIMENTS THAT HAVE MULTIDIR
count=1;
SubplotsPerFig=4;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

for i=1:NumBirds;
    
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        
        title([birdname ' - ' exptname]);
        xlabel('Days');
        ylabel('FF (hz minus baseline)');
        
        % === PLOT ACROSS DAYS LEARNING
        unique_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;  % all syls
        multidir_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls; % multidir syls
        
        % --- Plot all syls, weak color
        for j=1:length(unique_syls);
            syl= unique_syls{j};
            
            % skip if will plot as multidir syl below;
            if any(strcmp(multidir_syls, syl))==1;
                continue;
            end
            Y = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase; % day means
            Ysem = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).semFF;
            X = 1:length(Y); % days
            
            % similar or different?
            is_similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            
            
            if is_similar==1;
                %                 shadedErrorBar(X, Y, Ysem, {'Color', plotcols{j}},1);
                %                 hplot(j)=lt_plot(X, Y, {'Color', plotcols{j}});
                plot(X, Y, '-b', 'LineWidth',1.5);
            else
                plot(X, Y, '-r', 'LineWidth',1.5);
            end
            
            %             shadedErrorBar(X, Y, Ysem, {'Color', 'k'},1);
        end
        
        % --- Plot the multidir targets boldly
        plotcols=lt_make_plot_colors(length(multidir_syls), 0, 0);
        hplot=[];
        for j=1:length(multidir_syls);
            syl= multidir_syls{j};
            
            Y = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase; % day means
            Ysem = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).semFF;
            X = 1:length(Y); % days
            
            shadedErrorBar(X, Y, Ysem, {'Color', plotcols{j}},1);
            hplot(j)=lt_plot(X, Y, {'Color', plotcols{j}});
        end
        legend(hplot, multidir_syls);
        
        % ---- put lines for start and end of stuff
        % WN
        WNon_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        WNoff_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        line([WNon_ind-0.6 WNon_ind-0.6], ylim, 'Color','k', 'LineStyle' ,'--', 'LineWidth',2);
        line([WNoff_ind+0.6 WNoff_ind+0.6], ylim, 'Color','k', 'LineStyle' ,'--', 'LineWidth',2);
        
        % multidir
        MDonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
        MDoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
        line([MDonInd-0.4 MDonInd-0.4], ylim, 'Color','r', 'LineWidth',2);
        line([MDoffInd+0.4 MDoffInd+0.4], ylim, 'Color','r','LineWidth', 2);
        
        line(xlim, [0 0], 'LineStyle', '--');
        
    end
end

%% OVERLAY ALL TRAJECTORIES FROM ALL BIRDS (LOCKED TO START OF MULTIDIR LEARNING)

DayBinSize=PARAMS.global.MULTIDIR.DayBinSize;

plotcols=lt_make_plot_colors(TotalNumExperiments, 0, 0);
cc=1;

count=1;
SubplotsPerFig=12;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


for i=1:NumBirds;
    
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;        
        multidir_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls; % multidir syls
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname ' - ' exptname]);
        xlabel(['Days (first multidir WN day is day ' num2str(DayBinSize+1)])
        ylabel('FF (hz minus baseline)');
        % What days to plot? (start, N days before start multidir; end:
        % N days at end of multidir)
        DayBeforeMDon_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
        MDoff_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
        
        daystokeep=DayBeforeMDon_ind-DayBinSize+1:MDoff_ind; % will take N days BEFORE start. end at last day (not post).
        Ycomp={};
        for j=1:length(multidir_syls);
            syl=multidir_syls{j};
            
            % -- EXTRACT DATA
            Y = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase;
            Ysem= SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).semFF;
            
            % take multidir days
            Y=Y(daystokeep);
            Ysem=Ysem(daystokeep);
            
            % --- PLOT
            shadedErrorBar(1:length(Y), Y, Ysem, {'Color', plotcols{cc}},1);
            lt_plot(1:length(Y), Y, {'Color',plotcols{cc}});
            
            % --- SAVE DATA TO OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).AcrossEpoch.FFmean_MinusBase=Y;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).AcrossEpoch.FFsem_MinusBase=Ysem;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).AcrossEpoch.DayInds=daystokeep;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).AcrossEpoch.NumPreWNDays=DayBinSize;  
            
            % Save below to calculate separation across days
            Ycomp{j}=Y;
        end
        
        % ======== CALCULATE SEPARATION ACROSS DAYS
        if length(Ycomp)==2;
            SeparationAcrossDays_MinusBase=Ycomp{1}-Ycomp{2};
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSylData.AcrossEpoch.Separation_UsingHzMinusBaseline=SeparationAcrossDays_MinusBase;
        else
            % don't analyze separation for expts with >2 targets
        end
        % ==========================
        
        % lines
        line([DayBinSize+0.5 DayBinSize+0.5], ylim, 'Color','k');
        
        
        cc=cc+1; % for plot color
    end
end








%% PLOT (in terms of simply separation) (also get other measures of how different/similar they are, etc)

% plotcols=lt_make_plot_colors(TotalNumExperiments, 0, 0);
% cc=1;
% 
% count=1;
% SubplotsPerFig=12;
% subplotrows=4;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];


for i=1:NumBirds;
    
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        multidir_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls; % multidir syls
        
        
        % =========== What days to plot? (start, N days before start multidir; end:
        % N days at end of multidir)
        DayBeforeMDon_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
        MDoff_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
        
        DayBins_start=DayBeforeMDon_ind-DayBinSize+1:DayBeforeMDon_ind;
        DayBins_end=MDoff_ind-DayBinSize+1:MDoff_ind;
        
        DayWindowsList={DayBins_start, DayBins_end};
        
        % --- OUTPUT
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Multidir_DayBins_start=DayBins_start;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Multidir_DayBins_end=DayBins_end;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Multidir_DayWindowsList=DayWindowsList;
        
        % =========== COLLECT STATS
        if length(multidir_syls)==2;
            % 2 syls only
            
            for iii=1:length(DayWindowsList);
                
                dayInds=DayWindowsList{iii};
                
                ffmean_actual={};
                ffsem={};
                ffmean_minusbase={};
                for j=1:length(multidir_syls);
                    syl=multidir_syls{j};
                    
                    % -- COLLECT FFVALS (across days)
                    ffvals_actual=[];
                    ffvals_minus_base=[];
                    tvals=[];
                    for k=1:length(dayInds);
                        day=dayInds(k);
                        
                        ffvals_actual=[ffvals_actual cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{day})];
                        ffvals_minus_base=[ffvals_minus_base ffvals_actual-SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF];
                        
                        tvals=[tvals SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day}];
                    end
                    ffmean_actual{j}=mean(ffvals_actual);
                    ffsem{j}=lt_sem(ffvals_actual);
                    ffmean_minusbase{j}=mean(ffvals_minus_base);
                    % ------------------------
                    
                    % PUT INTO OUTPUT
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).DayWindow{iii}.ffvals_actual=ffvals_actual;
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).DayWindow{iii}.ffvals_minus_base=ffvals_minus_base;
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).DayWindow{iii}.tvals=tvals;
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).DayWindow{iii}.ffmean_actual=ffmean_actual{j};
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).DayWindow{iii}.ffsem=ffsem{j};
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(syl).DayWindow{iii}.ffmean_minusbase=ffmean_minusbase{j};
                end
                
                % =========== STATS RELATIVE BETWEEEN TWO SYLS
                Separation_UsingHzMinusBaseline=ffmean_minusbase{1}-ffmean_minusbase{2};
                Separation_UsingActualHz=ffmean_actual{1}-ffmean_actual{2};
                
                % --- PUT INTO OUTPUT
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSylData.DayWindow{iii}.Separation_UsingHzMinusBaseline=Separation_UsingHzMinusBaseline;
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSylData.DayWindow{iii}.Separation_UsingActualHz=Separation_UsingActualHz;
                
            end
            
            % ========== STATS RELATIVE (NOT RELATED TO SEPARATION)
            % --- Similarity (acoustic)
            fv1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(multidir_syls{1}).fv_baseline_zscore_mean;
            fv2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(multidir_syls{2}).fv_baseline_zscore_mean;
            
            acoustic_dist=sqrt(sum((fv1-fv2).^2));
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSyls_OtherStats.acoustic_dist=acoustic_dist;
            
            % --- Correlation
            corr=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(multidir_syls{1}).CORRELATIONS.song_by_song.corrcoeff_vs.(multidir_syls{2});
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSyls_OtherStats.correlation_SbyS=corr;
            
            % --- Magnitude of initial generalization
            FF=[];
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs~=2;            % not defined for experiments start started out targeting
                % both syllables
                FF(1)=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(multidir_syls{1}).DayWindow{1}.ffmean_minusbase;
                FF(2)=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(multidir_syls{2}).DayWindow{1}.ffmean_minusbase;
                
                % initial generalization
                [~, ind]=max(abs(FF)); % figure out which one was first target (it shows greater learning)
                
                first_targsyl_FF=FF(ind);
                other_syl_FF=FF(3-ind); % i.e. if first is 2, then this is 1.
                
                initial_generalization=other_syl_FF/first_targsyl_FF;
                
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSyls_OtherStats.initial_generalization=initial_generalization;
                
            % -------- How much does separation increase during bidir
            % learning?
            
            
            end
        else
            % more than 2 syls, analyze separately
        end
    end
end



% ================== PLOT (DOTS)
lt_figure; hold on;
title('Separation during bidir learning');

for i=1:NumBirds;
    
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        multidir_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls; % multidir syls
        
        if length(multidir_syls)~=2;
            
            % 2 syls only
            continue
        end
        
        X=[];
        Y=[];
        
        DayWindowsList=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.Multidir_DayWindowsList;
        for iii=1:length(DayWindowsList);
            
            dayInds=DayWindowsList{iii};
            
            % --- EXTRACT DATA
            Separation_UsingHzMinusBaseline=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSylData.DayWindow{iii}.Separation_UsingHzMinusBaseline;
            Separation_UsingActualHz=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSylData.DayWindow{iii}.Separation_UsingActualHz;
            
            X=[X iii];
            Y=[Y Separation_UsingHzMinusBaseline];
            
        end
        
        % If Y(end) is negative, then flip sign of all values in Y
        if Y(end)<0;
            Y=-1.*Y;
        end
        
        
        % === PLOT
        
        lt_plot(X, Y, {'LineStyle', '-'});
        
        % ==== PERFORM PERMUTATION TEST OF SIGNIFICANCE
        % method: [syl1_early syl1_late] [syl2_early syl2_late] - shuffle
        % each separately. then subtract (syl2 - syl1)ealry and
        % (syl2-syl1)late, then subtrtact those 2 things to get a
        % difference.
        % rewrite that at [x1 x2] [y1 y2];
        x1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.SylData.(multidir_syls{1}).DayWindow{1}.ffvals_minus_base;
        
        
        
        
    end
    
end

xlim([X(1)-1 X(2)+1]);


%% ================== PLOT SEPARATION ACROSS DAYS FOR ALL EXPERIMENTS ON ONE PLOT

count=1;
SubplotsPerFig=12;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:NumBirds;
    
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        multidir_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls; % multidir syls
        
        % figure stuff
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname ' - ' exptname]);
        xlabel(['Days (first multidir WN day is day ' num2str(DayBinSize+1) ')'])
        ylabel('Separation (hz)');
        
        if length(multidir_syls)>2;
            continue;
        end
        
        % ==== EXTRACT DATA
        Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSylData.AcrossEpoch.Separation_UsingHzMinusBaseline;
        
        if Y(end)<0;
            Y=-1.*Y;
        end
        
        % PLOT
        lt_plot(1:length(Y), Y, {'Color','b','LineStyle','-'});
        
        
            lt_plot_zeroline;        
            line([DayBinSize+0.5 DayBinSize+0.5], ylim, 'Color', 'k');
    
            
    end
end

lt_subtitle('Separation driven by 2-direction WN (using hz minus baseline)');


%% ================== PLOT SEPARATION ACROSS DAYS FOR ALL EXPERIMENTS ON ONE PLOT
% ORDER THE SUBPLOTS BASED ON SCORES LIKE CORRELATIONS (EASY TO DO WITH
% ACOUSTIC, ETC)

Y_all={};
Corr_all=[];
Birdnames={};
Exptnames={};
SylList={};
Acoustic_all=[];
InitialGeneralization=[];

count=1;
for i=1:NumBirds;
    
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        multidir_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls; % multidir syls
        
        if length(multidir_syls)>2;
            continue;
        end
        
        % ==== EXTRACT DATA
        Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSylData.AcrossEpoch.Separation_UsingHzMinusBaseline;
        
        Yall{count}=Y;
        
        % initial generalization
        target=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        syl_nottarg=multidir_syls{~strcmp(multidir_syls, target)};
        
        % --- OTHER THINGS
        Birdnames{count}=birdname;
        Exptnames{count}=exptname;
        SylList{count}=multidir_syls;
        
        
        InitialGeneralization(count)=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl_nottarg).LEARNING.consolid_start_rel_targ;
        
        
        
        % CORRELATIONS + ACOUSTIC
        corr = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSyls_OtherStats.correlation_SbyS;
        acoustic=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSyls_OtherStats.acoustic_dist;
        Corr_all(count)=corr;
        Acoustic_all(count)=acoustic;
        
        % --------
        count=count+1;
    end
end

% ====== PLOT (ordered by initial generalization)
[~, OrderOfPlot]=sort(Corr_all);
hplot=[];

count=1;
SubplotsPerFig=12;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(OrderOfPlot);
    IndToPlot=OrderOfPlot(i);
    
    
    % figure stuff
    [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
    %         title([birdname ' - ' exptname]);
    %         xlabel(['Days (first multidir WN day is day ' num2str(DayBinSize+1) ')'])
    %         ylabel('Separation (hz)');
   
    hplot(i)=gca;
    
    Y=Yall{IndToPlot};
    corr=Corr_all(IndToPlot);
    birdname=Birdnames{IndToPlot};
    exptname=Exptnames{IndToPlot};
    
    title([birdname '-' exptname '; corr: ' num2str(corr)])
    
    if Y(end)<0;
        Y=-1.*Y;
    end
    
    plot(1:length(Y), Y, '-ob');
end

linkaxes(hplot, 'y');

% ====== PLOT (ACOUSTIC)
[~, OrderOfPlot]=sort(Acoustic_all);
hplot=[];

count=1;
SubplotsPerFig=12;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(OrderOfPlot);
    IndToPlot=OrderOfPlot(i);
    
    
    [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
   
    hplot(i)=gca;
    
    Y=Yall{IndToPlot};
    acoustic=Acoustic_all(IndToPlot);
    birdname=Birdnames{IndToPlot};
    exptname=Exptnames{IndToPlot};
    
    title([birdname '-' exptname '; acou. dist.: ' num2str(acoustic)])
    
    if Y(end)<0;
        Y=-1.*Y;
    end
    
    plot(1:length(Y), Y, '-ob');
end

linkaxes(hplot, 'xy');
lt_subtitle('subplots ordered by acoustic distance');

% ====== PLOT (INIT GENERAL)
[~, OrderOfPlot]=sort(InitialGeneralization);
hplot=[];

count=1;
SubplotsPerFig=12;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:length(OrderOfPlot);
    IndToPlot=OrderOfPlot(i);
    
    
    [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
   
    hplot(i)=gca;
    
    Y=Yall{IndToPlot};
    acoustic=InitialGeneralization(IndToPlot);
    birdname=Birdnames{IndToPlot};
    exptname=Exptnames{IndToPlot};
    
    title([birdname '-' exptname '; init. general.: ' num2str(acoustic)])
    
    if Y(end)<0;
        Y=-1.*Y;
    end
    
    plot(1:length(Y), Y, '-ob');
end

linkaxes(hplot, 'x');
lt_subtitle('subplots ordered by initial genearlization');

%% ========= PLOT SEPARATION, ALL ON ONE PLOT
lt_figure; hold on;
title('Separation driven by bidir WN (using hz, minus baseline)');
xlabel(['Days (first multidir WN day is day ' num2str(DayBinSize+1) ')'])
ylabel('Separation (hz)');

plotcols=lt_make_plot_colors(TotalNumExperiments, 0, 0);
cc=1;

for i=1:NumBirds;
    
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        multidir_syls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls; % multidir syls
        
        if length(multidir_syls)>2;
            continue;
        end
        
        % ==== EXTRACT DATA
        Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_MultiDirLearning.ComparingSylData.AcrossEpoch.Separation_UsingHzMinusBaseline;
        X=1:length(Y);
        
        
        % remove nans
        inds=~isnan(Y);
        X=X(inds);
        Y=Y(inds);
       
        if Y(end)<0;
            Y=-1.*Y;
        end
        
        lt_plot(X, Y, {'Color',plotcols{cc},'LineStyle','-', 'LineWidth',2});
        
        cc=cc+1;
    end
end

lt_plot_zeroline;
line([DayBinSize+0.5 DayBinSize+0.5], ylim, 'Color', 'k');
