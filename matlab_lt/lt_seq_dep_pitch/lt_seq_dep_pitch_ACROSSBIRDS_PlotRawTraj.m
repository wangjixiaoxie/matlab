function lt_seq_dep_pitch_ACROSSBIRDS_PlotRawTraj(SeqDepPitch_AcrossBirds, PARAMS, BirdToPlot, ExptToPlot, SylsToPlot, use_rand_color,UseRecalcBaseline, ThreshChangeTime)

PlotLearnShift=1; % overlays text info
UseRecalcBaseline=1; % then limits to days used to calcualte deviation of learning

%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% Initiate plots
count=1;
SubplotsPerFig=4;
subplotrows=4;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];


if ~strcmp(BirdToPlot, 'all')
    % then user defined one experiment (figure out with one)
    
    for i=1:NumBirds;
        
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        if ~strcmp(birdname, BirdToPlot)
            continue;
        end
        
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            
            if ~strcmp(ExptToPlot, exptID);
                continue
            end
            
            % ===== THIS IS THE EXPERIMENT - PLOT
            if strcmp(SylsToPlot, 'all')
                SylsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            end
            
            
            plotcols=lt_make_plot_colors(length(SylsToPlot), 0,0);
            hplots=[];
            for j=1:length(SylsToPlot);
                syl=SylsToPlot{j};
                
                [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
                title(syl); xlabel('days'); ylabel('FF (hz)');
                
                NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals);
                
                
                % ----- plot color
                if use_rand_color==1;
                    PlotCol=plotcols{j};
                else
                    similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                    presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                    
                    % use functional scheme
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
                end
                
                hplots=[hplots hsplot];
                
                % ---- raw dat, all days
                ymeans=[];
                ystd=[];
                for day=1:NumDays;
                    if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day});
                        ymeans=[ymeans nan];
                        ystd=[ystd nan];
                        
                        continue
                    end
                    
                    firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                        tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
                        tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                        tvals=tvals.FinalValue;
                        ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{day};
                        
                    else
                        tvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day});
                        tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                        tvals=tvals.FinalValue;
                        ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{day});
                    end
                    
                    plot(tvals, ffvals, '.', 'Color',PlotCol);
                    
                    % collect
                    ymeans=[ymeans mean(ffvals)];
                    ystd=[ystd std(ffvals)];
                end
                
                % ----- overlay means
                lt_plot([1:NumDays]+0.7, ymeans, {'LineStyle','-','Color', PlotCol, 'Errors', ystd, 'MarkerSize',8});
                
                
                
                % ---- baseline mean
                if UseRecalcBaseline
                        BlineMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_RecalcDays;
                        BlineSTD=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_RecalcDays);
                else
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                        BlineMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                        BlineSTD=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                    else
                        BlineMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
                        BlineSTD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).stdFF;
                    end
                end
                
                shadedErrorBar(1:NumDays, ones(1, NumDays)*BlineMean, ones(1, NumDays)*BlineSTD, {'Color',PlotCol}, 1);
                
                % ---- WN ON/OFF
                WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                WNoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
                
                line([WNonInd WNonInd], ylim, 'Color','k', 'LineWidth',2)
                line([WNoffInd+1 WNoffInd+1], ylim, 'Color','k', 'LineWidth',2)
                
                % ---- OVERLAY WN THRESHOLDS IF DEFINED FOR THIS
                % EXPERIMENT\
                %
                %                 if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1
                ind1=find(strcmp(ThreshChangeTime, birdname));
                ind2=find(strcmp(ThreshChangeTime, exptID));
                ind3=find(strcmp(ThreshChangeTime, syl));
                
                if ~any(cellfun(@isempty, {ind1, ind2, ind3}))
                    indind=find(ind1+1 == ind2 & ind1+2 == ind3);
                    if ~isempty(indind)
                        threshchanges=ThreshChangeTime{ind1(indind)+3};
                    end
                    % plot those thresh changes
                    learndir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
                    firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                    for k=1:length(threshchanges)
                        
                        tmp=lt_convert_EventTimes_to_RelTimes(firstday, {threshchanges{k}{1}});
                        x1=tmp.FinalValue; % time of start
                        
                        if k<length(threshchanges)
                            tmp=lt_convert_EventTimes_to_RelTimes(firstday, {threshchanges{k+1}{1}});
                            x2=tmp.FinalValue; % time of start
                        else
                            Xlim=xlim;
                            x2=Xlim(2);
                            % then go to end of time
                        end
                        
                        thresh=threshchanges{k}{2};
                        line([x1 x2], [thresh thresh]);
                        
                    end
                end
                % ---------
                
                
                % ==== plot learning value
                learningZ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                lt_plot_annotation(1, ['shift (z score) = ' num2str(learningZ, '%3.2g')], PlotCol)
                
            end
            
            linkaxes(hplots, 'x')
            
        end
    end
else
    disp('PLOTTING all birds not coded yet!!!');
end



%% Initiate plots [CLEAN VERSION]
count=1;
SubplotsPerFig=4;
subplotrows=4;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];


if ~strcmp(BirdToPlot, 'all')
    % then user defined one experiment (figure out with one)
    
    for i=1:NumBirds;
        
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        if ~strcmp(birdname, BirdToPlot)
            continue;
        end
        
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptID=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            
            if ~strcmp(ExptToPlot, exptID);
                continue
            end
            
            % ===== THIS IS THE EXPERIMENT - PLOT
            if strcmp(SylsToPlot, 'all')
                SylsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            end
            
            
            plotcols=lt_make_plot_colors(length(SylsToPlot), 0,0);
            hplots=[];
            for j=1:length(SylsToPlot);
                syl=SylsToPlot{j};
                
                [fignums_alreadyused, hfigs, count, hsplot]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
                title(syl); xlabel('days'); ylabel('FF (hz)');
                hold on;
                
                NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals);
                
                
                % ----- plot color
                if use_rand_color==1;
                    PlotCol=plotcols{j};
                else
                    similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
                    presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                    
                    % use functional scheme
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
                end
                
                hplots=[hplots hsplot];
                
                % ---- raw dat, all days
                ymeans=[];
                ystd=[];
                for day=1:NumDays;
                    if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day});
                        ymeans=[ymeans nan];
                        ystd=[ystd nan];
                        
                        continue
                    end
                    
                    firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                        tvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{day};
                        tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                        tvals=tvals.FinalValue;
                        ffvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals_WithinTimeWindow{day};
                        
                    else
                        tvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).Tvals{day});
                        tvals=lt_convert_EventTimes_to_RelTimes(firstday, tvals);  % convert tvals to dayvals
                        tvals=tvals.FinalValue;
                        ffvals=cell2mat(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{day});
                    end
                    
                    plot(tvals, ffvals, '.', 'Color',PlotCol);
                    
                    % collect
                    ymeans=[ymeans mean(ffvals)];
                    ystd=[ystd std(ffvals)];
                end
                
                % ----- overlay means
                lt_plot([1:NumDays]+0.7, ymeans, {'LineStyle','-','Color', PlotCol, 'Errors', ystd, 'MarkerSize',8});
                
                
                
                % ---- baseline mean
                if UseRecalcBaseline==1
                    BlineMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_RecalcDays;
                    BlineSTD=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_RecalcDays);
                else
                    if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                        BlineMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                        BlineSTD=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                    else
                        BlineMean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF;
                        BlineSTD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).stdFF;
                    end
                end
                
                %                 shadedErrorBar(1:NumDays, ones(1, NumDays)*BlineMean, ones(1, NumDays)*BlineSTD, {'Color',PlotCol}, 1);
                line(xlim, [BlineMean BlineMean], 'Color','k', 'LineWidth',2);
                
                % ---- WN ON/OFF
                WNonInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
                WNoffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
                
                line([WNonInd WNonInd], ylim, 'Color','k', 'LineWidth',2)
                line([WNoffInd+1 WNoffInd+1], ylim, 'Color','k', 'LineWidth',2)
                
                
                % ==== plot learning value
                learningZ=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
                learningDay=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.dayIndsUsed;
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                    baseSTD=std(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                    
                else
                    baseSTD=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).stdFF;
                end
                learningHz=learningZ*baseSTD;
                Ylim=ylim;
                lt_plot_text(1, Ylim(2)-100, ['shift (z score) = ' num2str(learningZ, '%3.2g') '; ' num2str(learningHz, '%3.4g') 'hz'], PlotCol); hold on;
                
                
            end
            
            linkaxes(hplots, 'x')
            
        end
    end
else
    disp('PLOTTING all birds not coded yet!!!');
end


