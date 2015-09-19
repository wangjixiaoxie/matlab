function [PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANsequence(SeqDepPitch_AcrossBirds, PARAMS)
%% LT 9/1/15 - using regular expreissions, plots occurance rate of all regular expression motifs - i.e. sequence learning, and effect of musc


%% PARAMS

NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% FOR EACH EXPERIMENT, PLOT TIMELINE (ONE VAL PER DAY) (EACH EXPRESSION ONE DATAPOINT)
% DOES FOR BOTH LMAN PBS AND MUSC (AND FOR NON-LMAN BIRDS TOO)

for i=1:NumBirds;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        expression_names=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions;
        num_expressions=length(expression_names);
        
        
        % ========= FOR LMAN PBS (OR JUST NON-LMAN BIRDS)
        musc_field=''; % '' for PBS, '_MUSC' for musc
        day_data_field=['day_data' musc_field];
        NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.(day_data_field));
        
        NumRends_Expr_AcrossDays=nan(NumDays, num_expressions);
        
        for j=1:NumDays;
            % for each day, go through all the expressions and get %
            % for each
            
            % skip if no data today
            if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.(day_data_field){j});
                continue;
            end
            
            for k=1:num_expressions;
                
                num_rends=size(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.(day_data_field){j}.data_WithOutlier{k}.Final_ARRAYS.FFvals,1);
                NumRends_Expr_AcrossDays(j,k)=num_rends;
                
            end
            
            
            
        end
        
        % ==== STATS AND SAVE THE RENDS ACROSS DAYS
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob.expression_names=expression_names;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob.(day_data_field).rends=NumRends_Expr_AcrossDays;
        
        num_rends_tot=sum(NumRends_Expr_AcrossDays,2);
        num_rends_frac=NumRends_Expr_AcrossDays./repmat(num_rends_tot,1,num_expressions);
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob.(day_data_field).fraction=num_rends_frac;
        
        
        
        % ===================== DO THE SAME FOR LMAN INACTIVATION (IF EXISTS)
        musc_field='_MUSC'; % '' for PBS, '_MUSC' for musc
        day_data_field=['day_data' musc_field];
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr, day_data_field);
        
        NumRends_Expr_AcrossDays=nan(NumDays, num_expressions);
        
        for j=1:NumDays;
            % for each day, go through all the expressions and get %
            % for each
            
            % skip if no data today
            if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.(day_data_field){j});
                continue;
            end
            
            for k=1:num_expressions;
                
                num_rends=size(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.(day_data_field){j}.data_WithOutlier{k}.Final_ARRAYS.FFvals,1);
                NumRends_Expr_AcrossDays(j,k)=num_rends;
                
            end
            
        end
        
        % ==== STATS AND SAVE THE RENDS ACROSS DAYS
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob.expression_names=expression_names;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob.(day_data_field).rends=NumRends_Expr_AcrossDays;
        
        num_rends_tot=sum(NumRends_Expr_AcrossDays,2);
        num_rends_frac=NumRends_Expr_AcrossDays./repmat(num_rends_tot,1,num_expressions);
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob.(day_data_field).fraction=num_rends_frac;
            
        end 
    end
end




%% ===== PLOT (PBS AND MUSC)      

count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


for i=1:NumBirds;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        expression_names=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions;
        num_expressions=length(expression_names);
        NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.day_data);
                
        
        % ===== PLOT (PBS, OR NON-LMAN)
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([SeqDepPitch_AcrossBirds.birds{i}.birdname '-' SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID]);
        plotcols=lt_make_plot_colors(num_expressions,0,0);
        
        num_rends_frac=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob.day_data.fraction;
        hplot=[];
        for k=1:num_expressions;
            X=1:NumDays;
            Y=num_rends_frac(:,k);
            
            hplot(k)=lt_plot(X, Y, {'LineStyle','-','Color',plotcols{k}});
        end
        legend(hplot, expression_names);
        
        % ===== PLOT (MUSC, OVERLAYED)
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob, 'day_data_MUSC');
                num_rends_frac=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_Sequence.data.regular_expression_prob.day_data_MUSC.fraction;
        for k=1:num_expressions;
            X=1:NumDays;
            Y=num_rends_frac(:,k);
            
            plot(X, Y, 's-','Color',plotcols{k}, 'MarkerSize',9);
        end
        end
        
        
        % ===== LINES
        % ----------- Put lines for WN start and end, and baseline, and
        % other things
        try
            WNOnInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            WnOffInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
            line([WNOnInd-0.5 WNOnInd-0.5],ylim,'LineStyle','--','Color','r'); % minus 0.5 since the datapoint is on the day, so want the line to be before datapojnt.
            line([WnOffInd+0.5 WnOffInd+0.5],ylim,'LineStyle','--','Color','r')
            text(WNOnInd, 0, 'WN', 'Color' , 'r');
            
            baseline_last_day=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.BaselineDays(end);
            line([baseline_last_day+0.4 baseline_last_day+0.4],ylim,'LineStyle','--','Color','g')
            text(baseline_last_day, 5, 'Baseline end', 'Color' , 'g');
            
            % snapshot - important, snapshot{1} is end of single dir learning,
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0;
                
                try
                    snapshotfield=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
                    snapshotfield=snapshotfield{1};
                    snapshot_start=snapshotfield(5:13);
                    firstday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
                    snapshot_start_ind=lt_convert_EventTimes_to_RelTimes(firstday, {snapshot_start});
                    snapshot_start_ind=snapshot_start_ind.FinalValue;
                    snapshot_end_ind=snapshot_start_ind+2;
                    
                    line([snapshot_start_ind-0.4 snapshot_start_ind-0.4], ylim, 'Color', 'k');
                    line([snapshot_end_ind+0.4 snapshot_end_ind+0.4], ylim, 'Color', 'k');
                    text(snapshot_start_ind, 5, 'Snapshot', 'Color' , 'k');
                catch err
                    disp(['NO SNAPSHOT for bird ' num2str(i) ', expt ' num2str(ii) ', but should be']);
                end
            end
            
            % Put line for multidir learning if applicable
            try % FOR LMAN BIRDS - using notes in Plot Learning Params.
                
                date1_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds-1;
                date2_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_end_Inds;
                
                % --- plot
                line([date1_ind date1_ind], ylim, 'Color', 'k', 'LineStyle', '--');
                line([date2_ind date2_ind], ylim, 'Color', 'k', 'LineStyle', '--');
                text(date1_ind, 50, ['Day before start, and last day of, multidir, using syls: ' ...
                    [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls]], 'Color' , 'k');
                
            catch err
            end
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart');
                date1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart;
                tmp=lt_convert_EventTimes_to_RelTimes(FirstDate, {date1});
                date1_ind=tmp.JustDays_rel;
                
                date2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay;
                tmp=lt_convert_EventTimes_to_RelTimes(FirstDate, {date2});
                date2_ind=tmp.JustDays_rel;
                
                % --- output
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind=date1_ind;
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind=date2_ind;
                
                % --- plot
                line([date1_ind date1_ind], ylim, 'Color', 'k', 'LineStyle', '--');
                line([date2_ind date2_ind], ylim, 'Color', 'k', 'LineStyle', '--');
                text(date1_ind, 50, ['Day before start, and last day of, multidir, using syls: ' ...
                    [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls]], 'Color' , 'k');
            end
        catch err
            disp('NOTE: SKIPPED PLOTTING LINES AND SUCH');
        end
        
        % ---- what was target?
        target=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        lt_plot_text(0,0.9, ['target: ' target]);
        
        try 
            multi_targs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
        
        lt_plot_text(0,0.8, ['multidir targs: ' multi_targs]);
        catch err
        end
        
        
    end
end

lt_subtitle('Motif prob (circle,PBS; sq, MUSC)');


