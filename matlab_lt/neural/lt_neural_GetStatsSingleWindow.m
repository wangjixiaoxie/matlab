function OUTSTRUCT = lt_neural_GetStatsSingleWindow(Yspks, windowstart, windowend)

% Yspks; % cell, each array timpointw in sec of spikes
% windowstart; % in sec (relative to onset of data)
% windowend;

%% Collect numspikes acrioss trails
numtrials = length(Yspks);

NumSpksAll = [];
for i=1: numtrials
   
    numspks = sum(Yspks{i} > windowstart & Yspks{i} < windowend);
    
    NumSpksAll = [NumSpksAll numspks];
    
end

%% === output
OUTSTRUCT = struct;

OUTSTRUCT.NumSpksRAW = NumSpksAll;

OUTSTRUCT.NumSpksMean = mean(NumSpksAll);
OUTSTRUCT.NumSpksStd = std(NumSpksAll);
OUTSTRUCT.NumSpksSEM = lt_sem(NumSpksAll);

OUTSTRUCT.FrateMean = OUTSTRUCT.NumSpksMean * (1/(windowend - windowstart));
OUTSTRUCT.FrateStd = OUTSTRUCT.NumSpksStd * (1/(windowend - windowstart));
OUTSTRUCT.FrateSEM = OUTSTRUCT.NumSpksSEM * (1/(windowend - windowstart));
