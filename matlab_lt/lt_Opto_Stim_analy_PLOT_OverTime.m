function lt_Opto_Stim_analy_PLOT_OverTime(StatsStruct,FieldsToCheck,Params)
%% LT 1/13/15 -
% Inputs:
% StatsStruct - produced by lt_Opto_Stim_analy_PLOT_Compare2

%% PARAMS

NumFields=length(FieldsToCheck);
PlotColors=lt_make_plot_colors(NumFields,0);

TimeWindowFields=fieldnames(StatsStruct.(FieldsToCheck{1}).WINDOWED);
NumTimeWindows=length(TimeWindowFields);

SmthBin=20; % moving avg window

%% PITCH OVER TIME

% PITCH
figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    subplot(NumTimeWindows,1,i); hold on;
    title([timewindfield]);
    ylabel('Pitch (hz)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=FieldsToCheck{ii};
    
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.vals;
        
        plot(X,Y,'.','Color',PlotColors{ii});
        
        % plot running avg (bin 20)
        Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.vals,SmthBin);
        plot(X,Y,'-','Color',PlotColors{ii});
        
        
        % plot mean
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.mean;
        Ysem=StatsStruct.(fieldname).WINDOWED.(timewindfield).Pitch.SEM;
        h1(ii)=errorbar(X(end)+0.5,Y,Ysem,'o','Color',PlotColors{ii});
                
    end
end

legend(h1,FieldsToCheck)


% ENTROPY
figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    subplot(NumTimeWindows,1,i); hold on;
    title([timewindfield]);
    ylabel('Entropy (log, -inf to 0)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=FieldsToCheck{ii};
    
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.vals;
        
        plot(X,Y,'.','Color',PlotColors{ii});
        
        % plot running avg (bin 20)
        Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.vals,SmthBin);
        plot(X,Y,'-','Color',PlotColors{ii});
        
        
        % plot mean
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.mean;
        Ysem=StatsStruct.(fieldname).WINDOWED.(timewindfield).WEntropy.SEM;
        h1(ii)=errorbar(X(end)+0.5,Y,Ysem,'o','Color',PlotColors{ii});
                
    end
end

legend(h1,FieldsToCheck)


% AMPLITUDE
figure; hold on;
X=[];
for i=1:NumTimeWindows;
    timewindfield=TimeWindowFields{i};
    
    subplot(NumTimeWindows,1,i); hold on;
    title([timewindfield]);
    ylabel('Entropy (log, -inf to 0)');
    xlabel('Time (hrs)');
    
    for ii=1:NumFields;
        fieldname=FieldsToCheck{ii};
    
        
        X=StatsStruct.(fieldname).WINDOWED.(timewindfield).Time.hours;
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.vals_log;
        
        plot(X,Y,'.','Color',PlotColors{ii});
        
        % plot running avg (bin 20)
        Y=smooth(StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.vals_log,SmthBin);
        plot(X,Y,'-','Color',PlotColors{ii});
        
        
        % plot mean
        Y=StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.mean;
        Ysem=StatsStruct.(fieldname).WINDOWED.(timewindfield).Ampl.SEM;
        h1(ii)=errorbar(X(end)+0.5,Y,Ysem,'o','Color',PlotColors{ii});
                
    end
end

legend(h1,FieldsToCheck)


%% SAVE

try cd('OverTime');
catch err
    mkdir('OverTime');
    cd('OverTime');
end

savemultfigs



