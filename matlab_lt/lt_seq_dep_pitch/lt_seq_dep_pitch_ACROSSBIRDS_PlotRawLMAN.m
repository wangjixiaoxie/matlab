function lt_seq_dep_pitch_ACROSSBIRDS_PlotRawLMAN(SeqDepPitch_AcrossBirds, PARAMS, BirdToPlot, ExptToPlot, SylsToPlot, overlayMeans, plotRawFF, UseSylColors, flipsign, use_std, OverlayLMANStats, OverlayMUSC_days, plotLarge)

if ~exist('plotLarge','var');
    plotLarge=0
end

plotLMANalldata=1; % then not just in time window. also plots swiching times.

if plotRawFF
    flipsign=0; % since is raw, don't need to flip
end

%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% Initiate plots
SylsToPlot_orig=SylsToPlot;

% === if want to plot all birds, then loop thru all birds and expts
% ---- initual figure for this expt

for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    
    for ii=1:numexpts
        exptID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % --- plot this expt, unless BirdToPlot exists, in which case only
        % plot if matches target expt
        if ~isempty(BirdToPlot)
            % then only want to plot one expt, skip if this is not that
            % expt
            if ~strcmp(ExptToPlot, exptID) | ~strcmp(birdname, BirdToPlot)
                continue
            end
        else
            % skip this expt if is not LMAN
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0
                continue
            end
        end
        
        SylsToPlot=[];
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        if strcmp(SylsToPlot_orig, 'all')
            SylsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        elseif strcmp(SylsToPlot_orig, 'same');
            SylsToPlot=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS ...
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS];
            SylsToPlot=[targsyl SylsToPlot];
        else
            SylsToPlot=SylsToPlot_orig;
        end
        
        assert(~isempty(SylsToPlot), 'PROBLEM< PLOT WHICH SYLS?')
        
        
        
        % =============== PLOT FOR THIS EXPT
        count=1;
        subplotrows=2;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        if plotLarge==1
        subplotrows=2;
        subplotcols=1;
        end
            
        
        hplots=[];
        for j=1:length(SylsToPlot);
            syl=SylsToPlot{j};
            
            [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
            title([birdname '-' exptID '-' syl]); xlabel('days'); ylabel('FF (hz)');
            
            NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals);
            
            
            % ----- plot color
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            % use functional scheme
            if UseSylColors==1;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                    PlotCol='k';
                elseif similar & presim
                    PlotCol='b';
                elseif similar & ~presim
                    PlotCol='c';
                elseif ~similar & presim
                    PlotCol='r';
                elseif ~similar & ~presim
                    PlotCol='m';
                end
            else
                PlotCol='k';
            end
            
            hplots=[hplots hsplot];
            
            % ========== RAW DAT [PBS]
            ymeans=[];
            ysem=[];
            TvalEdgeAll=[];
            
            for day=1:NumDays;
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day});
                    ymeans=[ymeans nan];
                    ysem=[ysem nan];
                    TvalEdgeAll=[TvalEdgeAll nan];
                    continue
                end
                
                if     plotLMANalldata==1
                    % then get all data from this day, not just time window
                    tvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day});
                    ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{day});
                    basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
                else
                    % then only time window
                    tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
                    ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{day};
                    basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                end
                
                % -- raw or base subtracted?
                if plotRawFF==0
                    ffvals=ffvals-basemean;
                    basemean=0;
                end
                
                % --- convert tvals to day
                firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                tvals=tvals.FinalValue;
                
                % -- flip sign
                if flipsign==1;
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
                        ffvals=-1*ffvals;
                        basemean=basemean*-1;
                    end
                end
                
                plot(tvals, ffvals, 'o', 'Color',PlotCol);
                
                
                % collect
                ymeans=[ymeans mean(ffvals)];
                ysem=[ysem lt_sem(ffvals)];
                TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
                
            end
            
            % ----- overlay means
            if overlayMeans==1;
                lt_plot(TvalEdgeAll, ymeans, {'Errors',ysem, 'Color', PlotCol, 'MarkerSize',8});
            end
            % --- line for baseline
            line(xlim, [basemean basemean], 'Color',PlotCol);
            
            
            % =========== RAW DAT MUSC
            ymeans=[];
            ysem=[];
            TvalEdgeAll=[];
            
            for day=1:NumDays;
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day});
                    ymeans=[ymeans nan];
                    ysem=[ysem nan];
                    TvalEdgeAll=[TvalEdgeAll nan];
                    
                    continue
                end
                
                if     plotLMANalldata==1
                    % then get all data from this day, not just time window
                    tvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals{day});
                    ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals{day});
                    basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF;
                else
                    % then only time window
                    tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
                    ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
                    basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                end
                
                % - convert tvals to day
                tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                tvals=tvals.FinalValue;
                
                % -- raw or base subtracted?
                if plotRawFF==0
                    ffvals=ffvals-basemean;
                    basemean=0;
                end
                
                
                
                %                 tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
                %                 firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                %                 tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                %                 tvals=tvals.FinalValue;
                %
                %                 if plotRawFF==1
                %                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
                %                     basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                %
                %                 else
                %                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
                %                     basemean=0;
                %                 end
                
                if flipsign==1;
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
                        ffvals=-1*ffvals;
                        basemean=basemean*-1;
                    end
                end
                
                plot(tvals, ffvals, 's', 'Color','r');
                
                % ===== plot lines for MUSC switch + lag
                lagtime=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.Lag_time;
                muscSchedule=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds{day};
                
                muscStart=day+muscSchedule.start/24;
                muscEnd=day+muscSchedule.end/24;
                
                line([muscStart muscStart], ylim, 'Color','M');
                line([muscStart+lagtime/24 muscStart+lagtime/24], ylim, 'Color','R');
                line([muscEnd muscEnd], ylim, 'Color','r');
                
                % collect
                ymeans=[ymeans mean(ffvals)];
                ysem=[ysem lt_sem(ffvals)];
                TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
            end
            
            % ----- overlay means
            if overlayMeans==1;
                lt_plot(TvalEdgeAll, ymeans, {'Marker','s','Errors',ysem, 'Color', 'r', 'MarkerSize',8});
            end
            
            line(xlim, [basemean basemean], 'Color',PlotCol, 'LineStyle', '--');
            
            
            % ---- WN ON/OFF
            WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            WNoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            
            line([WNonInd WNonInd], ylim, 'Color','k', 'LineWidth',2)
            line([WNoffInd WNoffInd], ylim, 'Color','k', 'LineWidth',2)
            
            
            %                 % ==== plot learning value
            %                 learningZ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            %                 lt_plot_annotation(1, ['shift (z score) = ' num2str(learningZ, '%3.2g')], PlotCol)
            
            % =====
            line(xlim, [0 0], 'Color', PlotCol);
            
        end
        
        
        
        % ================================================== plot just means
        
        for j=1:length(SylsToPlot);
            syl=SylsToPlot{j};
            
            [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
            title([birdname '-' exptID '-'  syl]); xlabel('days'); ylabel('FF (hz)');
            
            NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals);
            
            
            % ----- plot color
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
            
            % use functional scheme
            if UseSylColors==1;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                    PlotCol='k';
                elseif similar & presim
                    PlotCol='b';
                elseif similar & ~presim
                    PlotCol='c';
                elseif ~similar & presim
                    PlotCol='r';
                elseif ~similar & ~presim
                    PlotCol='m';
                end
            else
                PlotCol='k';
            end
            
            hplots=[hplots hsplot];
            
            % ========== RAW DAT [PBS]
            ymeans=[];
            ysem=[];
            TvalEdgeAll=[];
            ystd=[];
            for day=1:NumDays;
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day});
                    ymeans=[ymeans nan];
                    ysem=[ysem nan];
                    TvalEdgeAll=[TvalEdgeAll nan];
                    ystd=[ystd nan];
                    
                    continue
                end
                
                tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
                firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                tvals=tvals.FinalValue;
                if plotRawFF==1
                    ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{day};
                    basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                else
                    ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
                    basemean=0;
                end
                
                if flipsign==1;
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
                        ffvals=-1*ffvals;
                        basemean=basemean*-1;
                    end
                end
                
                %                 plot(tvals, ffvals, 'o', 'Color',PlotCol);
                
                
                % collect
                ymeans=[ymeans mean(ffvals)];
                ysem=[ysem lt_sem(ffvals)];
                ystd=[ystd nanstd(ffvals)];
                TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
                
            end
            
            % ----- overlay means
            if use_std==1
                lt_plot(TvalEdgeAll, ymeans, {'Errors',ystd, 'Color', PlotCol, 'MarkerSize',8, 'LineStyle','-'});
            else
                lt_plot(TvalEdgeAll, ymeans, {'Errors',ysem, 'Color', PlotCol, 'MarkerSize',8, 'LineStyle','-'});
            end
            % --- line for baseline
            line(xlim, [basemean basemean], 'Color',PlotCol);
            
            
            % =========== RAW DAT MUSC
            ymeans=[];
            ysem=[];
            TvalEdgeAll=[];
            ystd=[];
            
            for day=1:NumDays;
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day});
                    ymeans=[ymeans nan];
                    ysem=[ysem nan];
                    ystd=[ystd nan];
                    TvalEdgeAll=[TvalEdgeAll nan];
                    
                    continue
                end
                
                tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
                firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                tvals=tvals.FinalValue;
                
                if plotRawFF==1
                    ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
                    basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                    
                else
                    ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
                    basemean=0;
                end
                
                if flipsign==1;
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
                        ffvals=-1*ffvals;
                        basemean=basemean*-1;
                    end
                end
                
                %                 plot(tvals, ffvals, 's', 'Color',PlotCol);
                
                % collect
                ymeans=[ymeans mean(ffvals)];
                ysem=[ysem lt_sem(ffvals)];
                ystd=[ystd std(ffvals)];
                TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
            end
            
            % ----- overlay means
            inds=~isnan(TvalEdgeAll);
            if use_std==1
                lt_plot(TvalEdgeAll(inds), ymeans(inds), {'LineStyle','--', 'Marker','s','Errors',ystd(inds), 'Color', 'r', 'MarkerSize',8});
            else
                lt_plot(TvalEdgeAll(inds), ymeans(inds), {'LineStyle','--', 'Marker','s','Errors',ysem(inds), 'Color', 'r', 'MarkerSize',8});
            end
            
            line(xlim, [basemean basemean], 'Color',PlotCol, 'LineStyle', '--');
            
            
            % ---- WN ON/OFF
            WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            WNoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            
            line([WNonInd WNonInd], ylim, 'Color','k', 'LineWidth',2)
            line([WNoffInd WNoffInd], ylim, 'Color','k', 'LineWidth',2)
            
            
            %                 % ==== plot learning value
            %                 learningZ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            %                 lt_plot_annotation(1, ['shift (z score) = ' num2str(learningZ, '%3.2g')], PlotCol)
            
            % =====
            line(xlim, [0 0], 'Color', PlotCol);
            
            
            % ==== OVERLAY LMAN LAERNING STATS?
            if OverlayLMANStats==1
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData, 'final_extracted_window')
                    
                    if isempty(OverlayMUSC_days)
                        % then use learning window already defined
                        ff_pbs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).meanFF_pbs;
                        ff_pbs_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).semFF_pbs;
                        ff_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).meanFF_musc;
                        ff_musc_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).semFF_musc;
                        
                    else
                        % use the user defined days
                        ffvals_pbs=[];
                        ffvals_musc=[];
                        for dayday=OverlayMUSC_days(1):OverlayMUSC_days(2)
                            
                            tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{dayday};
                            ffvals_pbs=[ffvals_pbs tmp];
                            
                            tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{dayday};
                            ffvals_musc=[ffvals_musc tmp];
                            
                        end
                        
                        ff_pbs=mean(ffvals_pbs);
                        ff_pbs_sem=lt_sem(ffvals_pbs);
                        ff_musc=mean(ffvals_musc);
                        ff_musc_sem=lt_sem(ffvals_musc);
                        
                    end
                    
                    if flipsign==1
                        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
                            ff_pbs=ff_pbs*-1;
                            ff_musc=ff_musc*-1;
                        end
                    end
                    
                    if plotRawFF==0
                        XLim=xlim;
                        shadedErrorBar(XLim, ones(1,2)*ff_pbs, ff_pbs_sem, {'LineStyle','--', 'Color','k'}, 1)
                        shadedErrorBar(XLim, ones(1,2)*ff_musc, ff_musc_sem, {'LineStyle','--', 'Color','r'}, 1)
                        %                     line(xlim, [ff_pbs ff_pbs], 'LineStyle','--', 'Color','k');
                        %                     line(xlim, [ff_pbs ff_pbs], 'LineStyle','--', 'Color','k');
                        %                     line(xlim, [ff_pbs ff_pbs], 'LineStyle','--', 'Color','k');
                        %                     line(xlim, [ff_musc ff_musc], 'LineStyle','--', 'Color','r');
                    end
                end
            end
            
        end
        
        linkaxes(hplots, 'x')
        
    end
end



%% =========== RUNNING CV if plotting raw LMAN (all datapoints even outside of time window), then plot running CV
if plotLMANalldata==1 & ~isempty(BirdToPlot)
    BinSize=15;
    % skip if trying to plot all birds, too many figures
    SylsToPlot_orig=SylsToPlot;
    
    % === if want to plot all birds, then loop thru all birds and expts
    % ---- initual figure for this expt
    
    for i=1:NumBirds
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        
        for ii=1:numexpts
            exptID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            % --- plot this expt, unless BirdToPlot exists, in which case only
            % plot if matches target expt
            if ~isempty(BirdToPlot)
                % then only want to plot one expt, skip if this is not that
                % expt
                if ~strcmp(ExptToPlot, exptID) | ~strcmp(birdname, BirdToPlot)
                    continue
                end
            else
                % skip this expt if is not LMAN
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0
                    continue
                end
            end
            
            SylsToPlot=[];
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            if strcmp(SylsToPlot_orig, 'all')
                SylsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            elseif strcmp(SylsToPlot_orig, 'same');
                SylsToPlot=[SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS ...
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS];
                SylsToPlot=[targsyl SylsToPlot];
            else
                SylsToPlot=SylsToPlot_orig;
            end
            
            assert(~isempty(SylsToPlot), 'PROBLEM< PLOT WHICH SYLS?')
            
            
            
            % =============== PLOT FOR THIS EXPT
            %         count=1;
            %         subplotrows=3;
            %         subplotcols=2;
            %         fignums_alreadyused=[];
            %         hfigs=[];
            %
            %
            %         hplots=[];
            for j=1:length(SylsToPlot);
                syl=SylsToPlot{j};
                
                [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
                title([birdname '-' exptID '-' syl]); xlabel('days'); ylabel('FF (hz)');
                
                NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals);
                
                
                % ----- plot color
                similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                
                % use functional scheme
                if UseSylColors==1;
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
                        PlotCol='k';
                    elseif similar & presim
                        PlotCol='b';
                    elseif similar & ~presim
                        PlotCol='c';
                    elseif ~similar & presim
                        PlotCol='r';
                    elseif ~similar & ~presim
                        PlotCol='m';
                    end
                else
                    PlotCol='k';
                end
                
                hplots=[hplots hsplot];
                
                % ========== RAW DAT [PBS]
                ymeans=[];
                ysem=[];
                TvalEdgeAll=[];
                
                for day=1:NumDays;
                    if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day});
                        ymeans=[ymeans nan];
                        ysem=[ysem nan];
                        TvalEdgeAll=[TvalEdgeAll nan];
                        continue
                    end
                    
                    if     plotLMANalldata==1
                        % then get all data from this day, not just time window
                        tvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day});
                        ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{day});
                        basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
                    else
                        % then only time window
                        tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
                        ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{day};
                        basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                    end
                    
                    % -- raw or base subtracted?
                    if plotRawFF==0
                        ffvals=ffvals-basemean;
                        basemean=0;
                    end
                    
                    % --- convert tvals to day
                    firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                    tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                    tvals=tvals.FinalValue;
                    
                    
                    % ------------ CONVERT TO RUNNING CV
                    tvalRun=lt_running_stats(tvals, BinSize);
                    ffvalRun=lt_running_stats(ffvals, BinSize);
                    
                    
                    
                    %                 % collect
                    %                 ymeans=[ymeans mean(ffvals)];
                    %                 ysem=[ysem lt_sem(ffvals)];
                    %                 TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
                    
                    
                    
                    % ===== IF TODAY HAS MUSC, COMBINE WITH PBS BEFORE PLOTTING.
                    % OTHERWISE PLOT JSUT PBS
                    
                    if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day});
                        
                        plot(tvalRun.Mean, ffvalRun.STD./ffvalRun.Mean, '.k');
                        
                    else
                        % --- CHANGE NAME OF PBS STUFF
                        tvalsPBS=tvals;
                        ffvalsPBS=ffvals;
                        
                        % --- combine data with MUSC
                        if     plotLMANalldata==1
                            % then get all data from this day, not just time window
                            tvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals{day});
                            ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals{day});
                            basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF;
                        else
                            % then only time window
                            tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
                            ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
                            basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                        end
                        
                        % - convert tvals to day
                        tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                        tvals=tvals.FinalValue;
                        
                        % --------------- COMBINE DATA
                        tvals=[tvals tvalsPBS];
                        ffvals=[ffvals ffvalsPBS];
                        
                        % ---- sort in temporal order
                        [tvals, inds]=sort(tvals);
                        ffvals=ffvals(inds);
                        
                        % ------------ CONVERT TO RUNNING CV
                        tvalRun=lt_running_stats(tvals, BinSize);
                        ffvalRun=lt_running_stats(ffvals, BinSize);
                        
                        plot(tvalRun.Mean, ffvalRun.STD./ffvalRun.Mean, 'k.');
                        
                        
                        % ===== plot lines for MUSC switch + lag
                        lagtime=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.Lag_time;
                        muscSchedule=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MuscimolSchedule_ByDayInds{day};
                        preSwitch=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.PBS_window(1);
                                                
                        muscStart=day+muscSchedule.start/24;
                        muscEnd=day+muscSchedule.end/24;
                        preSwitch=day+(muscSchedule.start/24)+(preSwitch/24);
                        
                        line([muscStart muscStart], ylim, 'Color','M');
                        line([muscStart+lagtime/24 muscStart+lagtime/24], ylim, 'Color','R');
                        line([muscEnd muscEnd], ylim, 'Color','r');
                        line([preSwitch preSwitch], ylim, 'Color','b');
                        
                    end
                    
                end
                
                
                %             % =========== RAW DAT MUSC
                %             ymeans=[];
                %             ysem=[];
                %             TvalEdgeAll=[];
                %
                %             for day=1:NumDays;
                %                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day});
                %                     ymeans=[ymeans nan];
                %                     ysem=[ysem nan];
                %                     TvalEdgeAll=[TvalEdgeAll nan];
                %
                %                     continue
                %                 end
                %
                %                 if     plotLMANalldata==1
                %                     % then get all data from this day, not just time window
                %                     tvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals{day});
                %                     ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals{day});
                %                     basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF;
                %                 else
                %                     % then only time window
                %                     tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
                %                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
                %                     basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                %                 end
                %
                %                 % - convert tvals to day
                %                 tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                %                 tvals=tvals.FinalValue;
                %
                % % -- raw or base subtracted?
                %                 if plotRawFF==0
                %                     ffvals=ffvals-basemean;
                %                     basemean=0;
                %                 end
                %
                %
                %
                % %                 tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
                % %                 firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                % %                 tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                % %                 tvals=tvals.FinalValue;
                % %
                % %                 if plotRawFF==1
                % %                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
                % %                     basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                % %
                % %                 else
                % %                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
                % %                     basemean=0;
                % %                 end
                %
                %                 if flipsign==1;
                %                     if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
                %                         ffvals=-1*ffvals;
                %                         basemean=basemean*-1;
                %                     end
                %                 end
                %
                %                 % ------------ CONVERT TO RUNNING CV
                %                 tvalRun=lt_running_stats(tvals, BinSize);
                %                 ffvalRun=lt_running_stats(ffvals, BinSize);
                %
                %                 plot(tvalRun.Mean, ffvalRun.STD./ffvalRun.Mean, 's', 'Color','r');
                %
                %
                %                 % collect
                % %                 ymeans=[ymeans mean(ffvals)];
                % %                 ysem=[ysem lt_sem(ffvals)];
                % %                 TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
                %             end
                %
                % %             % ----- overlay means
                % %             if overlayMeans==1;
                % %                 lt_plot(TvalEdgeAll, ymeans, {'Marker','s','Errors',ysem, 'Color', 'r', 'MarkerSize',8});
                % %             end
                % %
                % %             line(xlim, [basemean basemean], 'Color',PlotCol, 'LineStyle', '--');
                %
                %
                % ---- WN ON/OFF
                WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                WNoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
                
                line([WNonInd WNonInd], ylim, 'Color','k', 'LineWidth',2)
                line([WNoffInd WNoffInd], ylim, 'Color','k', 'LineWidth',2)
                
                
                %                 % ==== plot learning value
                %                 learningZ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                %                 lt_plot_annotation(1, ['shift (z score) = ' num2str(learningZ, '%3.2g')], PlotCol)
                
                % =====
                line(xlim, [0 0], 'Color', PlotCol);
                
            end
            
            
            
            % ================================================== plot just means
            
            
            linkaxes(hplots, 'x')
            
        end
    end
    
end

%% OLD VERSION, SAME AS ABOVE, EXCEPT HAS TO TAKE IN BIRD,. EXPT AS INPUT (WILL NOT GO THRU ALL)
% count=1;
% subplotrows=2;
% subplotcols=1;
% fignums_alreadyused=[];
% hfigs=[];
%
%
% for i=1:NumBirds;
%
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     if ~strcmp(birdname, BirdToPlot)
%         continue;
%     end
%
%     numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%
%     for ii=1:numexpts
%         exptID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%
%
%         if ~strcmp(ExptToPlot, exptID);
%             continue
%         end
%
%         % ===== THIS IS THE EXPERIMENT - PLOT
%         if strcmp(SylsToPlot, 'all')
%             SylsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         end
%
%
%         hplots=[];
%         for j=1:length(SylsToPlot);
%             syl=SylsToPlot{j};
%
%             [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
%             title(syl); xlabel('days'); ylabel('FF (hz)');
%
%             NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals);
%
%
%             % ----- plot color
%             similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
%             presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
%
%             % use functional scheme
%             if UseSylColors==1;
%                 if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
%                     PlotCol='k';
%                 elseif similar & presim
%                     PlotCol='b';
%                 elseif similar & ~presim
%                     PlotCol='c';
%                 elseif ~similar & presim
%                     PlotCol='r';
%                 elseif ~similar & ~presim
%                     PlotCol='m';
%                 end
%             else
%                 PlotCol='k';
%             end
%
%             hplots=[hplots hsplot];
%
%             % ========== RAW DAT [PBS]
%             ymeans=[];
%             ysem=[];
%             TvalEdgeAll=[];
%
%             for day=1:NumDays;
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day});
%                     ymeans=[ymeans nan];
%                     ysem=[ysem nan];
%                     TvalEdgeAll=[TvalEdgeAll nan];
%                     continue
%                 end
%
%                 tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
%                 firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
%                 tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
%                 tvals=tvals.FinalValue;
%                 if plotRawFF==1
%                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{day};
%                     basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
%                 else
%                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
%                     basemean=0;
%                 end
%
%                 if flipsign==1;
%                     if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
%                         ffvals=-1*ffvals;
%                         basemean=basemean*-1;
%                     end
%                 end
%
%                 plot(tvals, ffvals, 'o', 'Color',PlotCol);
%
%
%                 % collect
%                 ymeans=[ymeans mean(ffvals)];
%                 ysem=[ysem lt_sem(ffvals)];
%                 TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
%
%             end
%
%             % ----- overlay means
%             if overlayMeans==1;
%                 lt_plot(TvalEdgeAll, ymeans, {'Errors',ysem, 'Color', PlotCol, 'MarkerSize',8});
%             end
%             % --- line for baseline
%             line(xlim, [basemean basemean], 'Color',PlotCol);
%
%
%             % =========== RAW DAT MUSC
%             ymeans=[];
%             ysem=[];
%             TvalEdgeAll=[];
%
%             for day=1:NumDays;
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day});
%                     ymeans=[ymeans nan];
%                     ysem=[ysem nan];
%                     TvalEdgeAll=[TvalEdgeAll nan];
%
%                     continue
%                 end
%
%                 tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
%                 firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
%                 tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
%                 tvals=tvals.FinalValue;
%
%                 if plotRawFF==1
%                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
%                     basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
%
%                 else
%                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
%                     basemean=0;
%                 end
%
%                 if flipsign==1;
%                     if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
%                         ffvals=-1*ffvals;
%                         basemean=basemean*-1;
%                     end
%                 end
%
%                 plot(tvals, ffvals, 's', 'Color','r');
%
%                 % collect
%                 ymeans=[ymeans mean(ffvals)];
%                 ysem=[ysem lt_sem(ffvals)];
%                 TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
%             end
%
%             % ----- overlay means
%             if overlayMeans==1;
%                 lt_plot(TvalEdgeAll, ymeans, {'Marker','s','Errors',ysem, 'Color', 'r', 'MarkerSize',8});
%             end
%
%             line(xlim, [basemean basemean], 'Color',PlotCol, 'LineStyle', '--');
%
%
%             % ---- WN ON/OFF
%             WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
%             WNoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
%
%             line([WNonInd WNonInd], ylim, 'Color','k', 'LineWidth',2)
%             line([WNoffInd WNoffInd], ylim, 'Color','k', 'LineWidth',2)
%
%
%             %                 % ==== plot learning value
%             %                 learningZ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
%             %                 lt_plot_annotation(1, ['shift (z score) = ' num2str(learningZ, '%3.2g')], PlotCol)
%
%             % =====
%             line(xlim, [0 0], 'Color', PlotCol);
%
%         end
%
%         linkaxes(hplots, 'x')
%
%     end
% end
%
% %% plot just means
%
% for i=1:NumBirds;
%
%     birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
%     if ~strcmp(birdname, BirdToPlot)
%         continue;
%     end
%
%     numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
%
%     for ii=1:numexpts
%         exptID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
%
%
%         if ~strcmp(ExptToPlot, exptID);
%             continue
%         end
%
%         % ===== THIS IS THE EXPERIMENT - PLOT
%         if strcmp(SylsToPlot, 'all')
%             SylsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%         end
%
%
%         hplots=[];
%         for j=1:length(SylsToPlot);
%             syl=SylsToPlot{j};
%
%             [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
%             title(syl); xlabel('days'); ylabel('FF (hz)');
%
%             NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals);
%
%
%             % ----- plot color
%             similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
%             presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
%
%             % use functional scheme
%             if UseSylColors==1;
%                 if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target;
%                     PlotCol='k';
%                 elseif similar & presim
%                     PlotCol='b';
%                 elseif similar & ~presim
%                     PlotCol='c';
%                 elseif ~similar & presim
%                     PlotCol='r';
%                 elseif ~similar & ~presim
%                     PlotCol='m';
%                 end
%             else
%                 PlotCol='k';
%             end
%
%             hplots=[hplots hsplot];
%
%             % ========== RAW DAT [PBS]
%             ymeans=[];
%             ysem=[];
%             TvalEdgeAll=[];
%             ystd=[];
%             for day=1:NumDays;
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day});
%                     ymeans=[ymeans nan];
%                     ysem=[ysem nan];
%                     TvalEdgeAll=[TvalEdgeAll nan];
%                     ystd=[ystd nan];
%
%                     continue
%                 end
%
%                 tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
%                 firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
%                 tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
%                 tvals=tvals.FinalValue;
%                 if plotRawFF==1
%                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{day};
%                     basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
%                 else
%                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
%                     basemean=0;
%                 end
%
%                 if flipsign==1;
%                     if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
%                         ffvals=-1*ffvals;
%                         basemean=basemean*-1;
%                     end
%                 end
%
%                 %                 plot(tvals, ffvals, 'o', 'Color',PlotCol);
%
%
%                 % collect
%                 ymeans=[ymeans mean(ffvals)];
%                 ysem=[ysem lt_sem(ffvals)];
%                 ystd=[ystd nanstd(ffvals)];
%                 TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
%
%             end
%
%             % ----- overlay means
%             if use_std==1
%             lt_plot(TvalEdgeAll, ymeans, {'Errors',ystd, 'Color', PlotCol, 'MarkerSize',8, 'LineStyle','-'});
%             else
%             lt_plot(TvalEdgeAll, ymeans, {'Errors',ysem, 'Color', PlotCol, 'MarkerSize',8, 'LineStyle','-'});
%             end
%             % --- line for baseline
%             line(xlim, [basemean basemean], 'Color',PlotCol);
%
%
%             % =========== RAW DAT MUSC
%             ymeans=[];
%             ysem=[];
%             TvalEdgeAll=[];
%             ystd=[];
%
%             for day=1:NumDays;
%                 if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day});
%                     ymeans=[ymeans nan];
%                     ysem=[ysem nan];
%                     ystd=[ystd nan];
%                     TvalEdgeAll=[TvalEdgeAll nan];
%
%                     continue
%                 end
%
%                 tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{day};
%                 firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
%                 tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
%                 tvals=tvals.FinalValue;
%
%                 if plotRawFF==1
%                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
%                     basemean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
%
%                 else
%                     ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
%                     basemean=0;
%                 end
%
%                 if flipsign==1;
%                     if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
%                         ffvals=-1*ffvals;
%                         basemean=basemean*-1;
%                     end
%                 end
%
%                 %                 plot(tvals, ffvals, 's', 'Color',PlotCol);
%
%                 % collect
%                 ymeans=[ymeans mean(ffvals)];
%                 ysem=[ysem lt_sem(ffvals)];
%                 ystd=[ystd std(ffvals)];
%                 TvalEdgeAll=[TvalEdgeAll max(tvals)+0.02]; % Collect Edge of tvals for mean plot
%             end
%
%             % ----- overlay means
%             inds=~isnan(TvalEdgeAll);
%             if use_std==1
%             lt_plot(TvalEdgeAll(inds), ymeans(inds), {'LineStyle','--', 'Marker','s','Errors',ystd(inds), 'Color', 'r', 'MarkerSize',8});
%             else
%             lt_plot(TvalEdgeAll(inds), ymeans(inds), {'LineStyle','--', 'Marker','s','Errors',ysem(inds), 'Color', 'r', 'MarkerSize',8});
%             end
%
%             line(xlim, [basemean basemean], 'Color',PlotCol, 'LineStyle', '--');
%
%
%             % ---- WN ON/OFF
%             WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
%             WNoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
%
%             line([WNonInd WNonInd], ylim, 'Color','k', 'LineWidth',2)
%             line([WNoffInd WNoffInd], ylim, 'Color','k', 'LineWidth',2)
%
%
%             %                 % ==== plot learning value
%             %                 learningZ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
%             %                 lt_plot_annotation(1, ['shift (z score) = ' num2str(learningZ, '%3.2g')], PlotCol)
%
%             % =====
%             line(xlim, [0 0], 'Color', PlotCol);
%
%
%             % ==== OVERLAY LMAN LAERNING STATS?
%             if OverlayLMANStats==1
%                 if isempty(OverlayMUSC_days)
%                     % then use learning window already defined
%                     ff_pbs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).meanFF_pbs;
%                     ff_pbs_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).semFF_pbs;
%                     ff_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).meanFF_musc;
%                     ff_musc_sem=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(syl).semFF_musc;
%
%                 else
%                     % use the user defined days
%                     ffvals_pbs=[];
%                     ffvals_musc=[];
%                     for dayday=OverlayMUSC_days(1):OverlayMUSC_days(2)
%
%                         tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_DevFromBase_WithinTimeWindow{dayday};
%                         ffvals_pbs=[ffvals_pbs tmp];
%
%                         tmp=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{dayday};
%                         ffvals_musc=[ffvals_musc tmp];
%
%                     end
%
%                     ff_pbs=mean(ffvals_pbs);
%                     ff_pbs_sem=lt_sem(ffvals_pbs);
%                     ff_musc=mean(ffvals_musc);
%                     ff_musc_sem=lt_sem(ffvals_musc);
%
%                 end
%
%                 if flipsign==1
%                     if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir==-1;
%                         ff_pbs=ff_pbs*-1;
%                         ff_musc=ff_musc*-1;
%                     end
%                 end
%
%                 if plotRawFF==0
%                     XLim=xlim;
%                     shadedErrorBar(XLim, ones(1,2)*ff_pbs, ff_pbs_sem, {'LineStyle','--', 'Color','k'}, 1)
%                     shadedErrorBar(XLim, ones(1,2)*ff_musc, ff_musc_sem, {'LineStyle','--', 'Color','r'}, 1)
% %                     line(xlim, [ff_pbs ff_pbs], 'LineStyle','--', 'Color','k');
% %                     line(xlim, [ff_pbs ff_pbs], 'LineStyle','--', 'Color','k');
% %                     line(xlim, [ff_pbs ff_pbs], 'LineStyle','--', 'Color','k');
% %                     line(xlim, [ff_musc ff_musc], 'LineStyle','--', 'Color','r');
%                 end
%
%             end
%
%         end
%
%         linkaxes(hplots, 'x')
%
%     end
% end
%
%

