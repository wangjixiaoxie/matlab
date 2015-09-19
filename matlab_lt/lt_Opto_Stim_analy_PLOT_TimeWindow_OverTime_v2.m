function [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_OverTime_v2(StatsStruct, Params,saveON,KeepOutliers)
%% LT 4/7/15 - v2 - overwrites structures to, save. changes dir format. better representation of running avg.
% KeepOutliers=0 (plots without Outliers), 1( plots with).

%% LT 1/13/15 -
% Inputs:
% StatsStruct, Params - produced by lt_Opto_Stim_analy_PLOT_TimeWindow
% saveON=1, will save, if 0, not
%% Auto PARAMS

NumFields=length(Params.FieldsToCheck);
PlotColors=lt_make_plot_colors(NumFields,0);

NumTimeWindows=length(Params.TimeWindowList);

%% PITCH OVER TIME
TimeWindowFields=Params.TimeField;

% determine subplot sizes
[~, row_plots, col_plots]=lt_get_subplot_size(NumTimeWindows,NumTimeWindows);


% ======================== PITCH
lt_figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    lt_subplot(row_plots,col_plots,i); hold on;
    title([timewindfield]);
    ylabel('Pitch (hz)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=Params.FieldsToCheck{ii};
        
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.vals;
        
        % outliers
        if KeepOutliers==0;
           X(StatsStruct.(fieldname).OutlierInds_AllWindsCombn)=[];
           Y(StatsStruct.(fieldname).OutlierInds_AllWindsCombn)=[];           
        end
        
        lt_plot(X,Y,{'MarkerSize',4,'Color',PlotColors{ii}});
        
        % sort to plot in order
        [X,inds]=sort(X);
        Y=Y(inds);
        
        % plot running avg (bin 20)
        if length(X)>Params.SmthBin;% only do if have enough data
            Y_sm=lt_running_stats(Y,Params.SmthBin);
            X_sm=lt_running_stats(X,Params.SmthBin);
            
            shadedErrorBar(X_sm.Median,Y_sm.Mean,Y_sm.SEM,{'-','LineWidth',2,'Color',PlotColors{ii}},1);
            %         Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.vals,Params.SmthBin);
            %         plot(X,Y,'-','Color',PlotColors{ii});
        end
        
        % plot mean
        Ymean=mean(Y);
        Ysem=lt_sem(Y);
        h1(ii)=errorbar(X(end)+0.1,Ymean,Ysem,'o','Color',PlotColors{ii});
        
    end
end

legend(h1,Params.FieldsToCheck)
% linkaxes(h1,'x')

% ================================== ENTROPY
lt_figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    lt_subplot(row_plots,col_plots,i); hold on;

    title([timewindfield]);
    ylabel('Entropy (log, -inf to 0)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=Params.FieldsToCheck{ii};
        
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.vals;
               
        % outliers
        if KeepOutliers==0;
           X(StatsStruct.(fieldname).OutlierInds_AllWindsCombn)=[];
           Y(StatsStruct.(fieldname).OutlierInds_AllWindsCombn)=[];           
        end

        lt_plot(X,Y,{'MarkerSize',4,'Color',PlotColors{ii}});
        
        % sort to plot in order
        [X,inds]=sort(X);
        Y=Y(inds);

        %         % plot running avg (bin 20)
        %         Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.vals,Params.SmthBin);
        %         plot(X,Y,'-','Color',PlotColors{ii});
        
        if length(X)>Params.SmthBin;% only do if have enough data
            Y_sm=lt_running_stats(Y,Params.SmthBin);
            X_sm=lt_running_stats(X,Params.SmthBin);
            
            shadedErrorBar(X_sm.Median,Y_sm.Mean,Y_sm.SEM,{'-','LineWidth',2,'Color',PlotColors{ii}},1);
            %         Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.vals,Params.SmthBin);
            %         plot(X,Y,'-','Color',PlotColors{ii});
        end
        
        
        % plot mean
        Ymean=mean(Y);
        Ysem=lt_sem(Y);
        h1(ii)=errorbar(X(end)+0.1,Ymean,Ysem,'o','Color',PlotColors{ii});
        
    end
end

legend(h1,Params.FieldsToCheck)
% linkaxes(h1,'x')


% ========================================= AMPLITUDE
lt_figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    lt_subplot(row_plots,col_plots,i); hold on;

    title([timewindfield]);
    ylabel('Amplitude (log, arbitrary)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=Params.FieldsToCheck{ii};
        
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.vals_log;
        
        % outliers
        if KeepOutliers==0;
           X(StatsStruct.(fieldname).OutlierInds_AllWindsCombn)=[];
           Y(StatsStruct.(fieldname).OutlierInds_AllWindsCombn)=[];           
        end

        lt_plot(X,Y,{'MarkerSize',4,'Color',PlotColors{ii}});
        
        % sort to plot in order
        [X,inds]=sort(X);
        Y=Y(inds);

        % plot running avg (bin 20)
        %         Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.vals_log,Params.SmthBin);
        %         plot(X,Y,'-','Color',PlotColors{ii});
        
        if length(X)>Params.SmthBin;% only do if have enough data
            Y_sm=lt_running_stats(Y,Params.SmthBin);
            X_sm=lt_running_stats(X,Params.SmthBin);
            
            shadedErrorBar(X_sm.Median,Y_sm.Mean,Y_sm.SEM,{'-','LineWidth',2,'Color',PlotColors{ii}},1);
            %         Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.vals,Params.SmthBin);
            %         plot(X,Y,'-','Color',PlotColors{ii});
        end
        
        
        % plot mean
        Ymean=mean(Y);
        Ysem=lt_sem(Y);
        h1(ii)=errorbar(X(end)+0.1,Ymean,Ysem,'o','Color',PlotColors{ii});
        
    end
end

legend(h1,Params.FieldsToCheck)
% linkaxes(h1,'x')


%% SAVE

if saveON==1;
    disp('Saving...')
    
    tstamp=lt_get_timestamp(0);
    cd(Params.savefolder);
    
    % save
    save('StatsStruct','StatsStruct');
    save('Params','Params');
    
    % make note
    DoneNote=['DONE_OverTime ' tstamp '.txt'];
    fid1=fopen(DoneNote,'w');
    fclose(fid1);
    
    % save figs
    try
        cd('FIGURES/OverTime');
    catch err
        mkdir('FIGURES/OverTime');
        cd('FIGURES/OverTime');
    end
    
    lt_save_all_figs;
    cd('../../');
    
end

disp('Done!');

