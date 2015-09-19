% LT 1/17/14 - use to load a structure compiled using
% lt_all_days_various_calculations program, which compiles analyses over
% days.

function [all_days_various]=lt_load_all_days_various_structure(bluejay_num, birdname, date_time)

% bluejay_num - e.g. bluejay2 would be 2.
% birdname - e.g. 'bk51bk59'
% date_time (e.g. '17Jan2014_1545'), that is identifier of the structure.

curr_dir=pwd;
cd(['/bluejay' num2str(bluejay_num) '/lucas/birds/' birdname '/all_days_various_calculations/'])
filename=dir(['all_days_various_' date_time '*']);
load(filename.name)
cd(curr_dir)

end
