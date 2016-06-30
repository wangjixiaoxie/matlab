function lt_seq_dep_pitch_ACROSSBIRDS_PlotAllTraj(SeqDepPitch_AcrossBirds, PARAMS, Absolute_FF)
%% LT 10/2/15 - Plots all learning trajectories for experiments in the SeqDepPitch_AcrossBirds structure
% Absolute_FF=1; % plots absolute FF

if ~exist('Absolute_FF', 'var');
    Absolute_FF=0;
end

%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% Plot learning across days - zscored

count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        plotcols=lt_make_plot_colors(length(syls_unique), 0, 0);
        hplot=[];
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            
            % ================ GATHER Z-SCORE LEARNING
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
                Zvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore_WithinTimeWindow;
                Y=Zvals;
                % ===== GATHER ABSOLUTE FF IF REQUIRED
                if Absolute_FF==1;
                    ffvals= SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_DevFromBase_WithinTimeWindow;
                    Y=ffvals;
                end
                
            else
                Zvals=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF_zscore;
                Y=Zvals;
                % ===== GATHER ABSOLUTE FF IF REQUIRED
                if Absolute_FF==1;
                    ffvals= SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).meanFF;
                    Y=ffvals;
                end
            end
            
            % +++++++++++++++++++
            % ==== PLOT
            plot(1:length(Y), Y, 'Color',plotcols{j});
            hplot(j)=lt_plot(1:length(Y), Y, {'Color',plotcols{j}});
        end
        
        legend(hplot, syls_unique);
        line(xlim, [1.5 1.5])
        line(xlim, [-1.5 -1.5])
        
        line(xlim, [2 2])
        line(xlim, [-2 -2])
        
        
        % ============== annotate with lines for stuff
        Ylim=ylim;
        Yrange=Ylim(2)-Ylim(1);
        Xlim=xlim;
        
        % WN on and off
        WNonday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        line([WNonday-0.5 WNonday-0.5], ylim, 'Color' ,'r');
        WNoffday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        line([WNoffday+0.5 WNoffday+0.5],ylim,'Color','r');
        lt_plot_text(WNonday, Ylim(1), 'WNon/off','r')
        
        % day I called consol start
        try
            day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
            day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolEndInd;
            
            line([day1-0.5 day1-0.5], ylim);
            line([day2+0.5 day2+0.5], ylim);
            
            lt_plot_text(day1, Ylim(1), 'consol','b')
        catch err
        end
        
        % multidir
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart');
            date1_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
            date2_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
            
            line([date1_ind+0.4 date1_ind+0.4], ylim, 'Color', 'k');
            line([date2_ind+0.4 date2_ind+0.4], ylim, 'Color', 'k');
            lt_plot_text(date1_ind, Ylim(1), 'multidir','k')
        end
        
        % LMAN inactivation days
        if  SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            muscdays_inds=find(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow));
            
            % plot
            plot(muscdays_inds, Ylim(1)+Yrange/6, '^k', 'MarkerSize', 8);
            plot(muscdays_inds, Ylim(2)-Yrange/6, 'vk', 'MarkerSize', 8);
            
            %             for k=1:length(muscdays_inds);
            % %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'Color',[0.7 0.7 0.7]);
            %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'LineStyle','--','Color',[0.3 0.3 0.3]);
            %             end
            
            lt_plot_text(1, Ylim(1)+Yrange/8, 'arrowhead: MUSC', 'k');
            
                        
        end
        
        % annotate target learning direction (as detected automatically)
        targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        lt_plot_text(2, Ylim(2)-Yrange/8, ['targ learn dir: ' num2str(targ_learn_dir)], 'r');
        
    end
end
















