function all_days_all_analysis=lt_load_all_days_all_analysis(ADAAindex)
% e.g. all_days_all_analysis=lt_load_all_days_all_analysis('17May2014_2020');

% RUN this from the bird directory
cd all_days_all_analysis/
fname=dir(['all_days_all_analysis_' ADAAindex '*']);
load(fname.name);
cd ..
