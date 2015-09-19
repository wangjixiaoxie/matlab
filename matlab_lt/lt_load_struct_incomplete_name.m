%% LT 6/2/14 - load a structure if you know the beginning of the name (but end variable), plus an asterisk
% e.g. [transitions_all_days, fname_complete]=lt_load_struct_incomplete_name('/bluejay2/lucas/birds/gr66gr43/all_days_transitions_ContextSeq2_analysis/transitions_all_days_01Jun2014_to_01Jun2014_div_a*')

function [transitions_all_days, fname_complete]=lt_load_struct_incomplete_name(incomplete_fname_asterisk);

fname_complete=ls(incomplete_fname_asterisk);
load(fname_complete(1:end-1))

% Alternatively do this, but it has to be run in the corrext dir.
% fname_complete=fname.name;
% load(fname_complete);

end
