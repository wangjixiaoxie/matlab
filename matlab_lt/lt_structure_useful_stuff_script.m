%% have array of structures, want to add one field without looping thru array inds
% if ss is like ss(1).field..., then do following [if value9 is same for
% all entries
[ss.field9]=deal( 'value9')

% if want to assign different value to each entry (e.g. ss(1).test=434;
% ss(2).test=545, do following
% if FFvals is array of numbers
c=num2cell(FFvals);
[data.test6]=deal(c{:});
% gives you data(1).test6 =FFvals(1), and so on...




