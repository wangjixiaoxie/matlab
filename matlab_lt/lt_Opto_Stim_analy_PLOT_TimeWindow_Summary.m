%% LT 2/16/15 - takes output of lt_Opto_Stim_analy_PLOT_TimeWindow and displays summary stats

function lt_Opto_Stim_analy_PLOT_TimeWindow_Summary(StatsStruct,Params);

% GET PARAMS
NumFields=length(Params.FieldsToCheck); % which trial classes? (e.g. stim vs. no stim);
TimeWindFields=Params.TimeField;

% RUN
disp(' ');
disp('DISPLAYING SUMMARY STATS:');

for i=1:NumFields; % all trial types
    fieldname=Params.FieldsToCheck{i};
    disp(' ');
    disp(['TRIAL TYPE: ' fieldname]);
    
    for ii=1:length(TimeWindFields);
        timefield=TimeWindFields{ii};
        disp(['   TIME BIN: ' timefield]);
        disp(['      Pitch (mean)= ' num2str(StatsStruct.(fieldname).WINDOWED.(timefield).Pitch.mean)]);
        disp(['      Pitch (SEM)= ' num2str(StatsStruct.(fieldname).WINDOWED.(timefield).Pitch.SEM)]);
        disp(['      Pitch (SD)= ' num2str(StatsStruct.(fieldname).WINDOWED.(timefield).Pitch.SD)]);
        disp(['      N= ' num2str(StatsStruct.(fieldname).WINDOWED.(timefield).n)]);
        
        try % because might not have median field
            disp(['      Pitch (median)= ' num2str(StatsStruct.(fieldname).WINDOWED.(timefield).Pitch.median)]);
        catch err
        end
        
        
    end
    
end

