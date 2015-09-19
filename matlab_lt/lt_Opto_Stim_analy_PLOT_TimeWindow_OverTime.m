function [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_OverTime(StatsStruct, Params,saveON)
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


% PITCH
figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    subplot(row_plots,col_plots,i); hold on;
    title([timewindfield]);
    ylabel('Pitch (hz)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=Params.FieldsToCheck{ii};
    
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.vals;
        
        plot(X,Y,'.','Color',PlotColors{ii});
        
        % plot running avg (bin 20)
        Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.vals,Params.SmthBin);
        plot(X,Y,'-','Color',PlotColors{ii});
        
        
        % plot mean
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.mean;
        Ysem=StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.SEM;
        h1(ii)=errorbar(X(end)+0.1,Y,Ysem,'o','Color',PlotColors{ii});
                
    end
end

legend(h1,Params.FieldsToCheck)
% linkaxes(h1,'x')

% ENTROPY
figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    subplot(row_plots,col_plots,i); hold on;
    title([timewindfield]);
    ylabel('Entropy (log, -inf to 0)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=Params.FieldsToCheck{ii};
    
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.vals;
        
        plot(X,Y,'.','Color',PlotColors{ii});
        
        % plot running avg (bin 20)
        Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.vals,Params.SmthBin);
        plot(X,Y,'-','Color',PlotColors{ii});
        
        
        % plot mean
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.mean;
        Ysem=StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.SEM;
        h1(ii)=errorbar(X(end)+0.1,Y,Ysem,'o','Color',PlotColors{ii});
                
    end
end

legend(h1,Params.FieldsToCheck)
% linkaxes(h1,'x')


% AMPLITUDE
figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    subplot(row_plots,col_plots,i); hold on;
    title([timewindfield]);
    ylabel('Amplitude (log, arbitrary)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=Params.FieldsToCheck{ii};
    
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.vals_log;
        
        plot(X,Y,'.','Color',PlotColors{ii});
        
        % plot running avg (bin 20)
        Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.vals_log,Params.SmthBin);
        plot(X,Y,'-','Color',PlotColors{ii});
        
        
        % plot mean
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.mean;
        Ysem=StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.SEM;
        h1(ii)=errorbar(X(end)+0.1,Y,Ysem,'o','Color',PlotColors{ii});
                
    end
end

legend(h1,Params.FieldsToCheck)
% linkaxes(h1,'x')


%% SAVE
if saveON==1;
tstamp=lt_get_timestamp(0);
sdir=[Params.PLOT_TimeWindow_savedir '/OverTime_' tstamp];
mkdir(sdir);
cd(sdir);

Params.PLOT_TimeWindow_OverTime_savedir=sdir;

savemultfigs
save('StatsStruct','StatsStruct');
save('Params','Params');
end


