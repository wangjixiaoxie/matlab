function lt_neural_v2_ListAllNeurons
%% LT 2/6/17 - lists all neurons that have been saved to summarystruct

% % load
% SummaryStruct = lt_neural_v2_LoadSummary;
% 
% ==== LIST
[NeuronDatabase] = lt_neural_v2_ConvertSummary2Database(''); % convert to database (then can make table)

% TABLE
NeuronDatabase_table=struct2table(NeuronDatabase.neurons);
summary(NeuronDatabase_table);
disp(NeuronDatabase_table)

f = figure('Position',[200 200 1000 500]);
dat=table2cell(NeuronDatabase_table);
cnames=fieldnames(NeuronDatabase_table);
t = uitable('Parent', f, 'Data',dat,'ColumnName',cnames, 'Position',[0 0 900 450]);

% === if there are any empty neurons, then move another neuron
