function SummaryStruct = lt_neural_v2_LoadSummary()

%% lt 2/5/17 - load summary struct. doesn't need inputs.

targdir = ['/bluejay5/lucas/analyses/neural/'];
clear SummaryStruct
load([targdir 'SummaryStruct']);
